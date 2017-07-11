//Added by Karel Adamek
//#define MSD_DEBUG

#include "headers/params.h"
#include "device_MSD_limited_kernel.cu"

int Choose_x_dim(int grid_dim){
	int seive[15] =	{ 32, 31, 29, 23, 19, 17, 16, 13, 11, 8, 7, 5, 4, 3, 2 };

	int f, nRest, nBlocks, N, N_accepted;

	N = 1;
	N_accepted = 1;
	for (int i = 0; i < 4; i++)	{
		for (f = 0; f < 15; f++) {
			nBlocks = grid_dim / seive[f];
			nRest = grid_dim - nBlocks*seive[f];
			if (nRest == 0) {
				N_accepted = N_accepted*N;
				N = seive[f];
				break;
			}
		}
		if (( N_accepted*N ) > 32 || N == 1)
			return ( N_accepted );
		grid_dim = grid_dim / N;
	}

	return ( N_accepted );
}

int Choose_y_dim(int grid_dim){
	int seive[5] = { 32, 16, 8, 4, 2 };

	int f, nRest, nBlocks, N;

	N = 1;
	for (f = 0; f < 5; f++) {
		nBlocks = grid_dim / seive[f];
		nRest = grid_dim - nBlocks*seive[f];
		if (nRest == 0) {
			N = seive[f];
			break;
		}
	}

	return ( N );
}

/*
int Choose_divider(int number, int max_divider){
	int seive[12]={2, 3, 4, 5, 7, 11, 13, 17, 19, 23, 29, 31};
	int f, nRest, nBlocks, N, N_accepted;
	
	N=1;N_accepted=1;
	do {
		N=1;
		for(f=0; f<12; f++){
			nBlocks=number/seive[f];
			nRest=number - nBlocks*seive[f];
			if(nRest==0) {
				N=seive[f];
				N_accepted=N_accepted*N;
				break;
			}
		}
		number=number/N;
	} while ( (N_accepted)<=max_divider && N>1 );
	
	return(N_accepted/N);
}
*/

void MSD_limited_init() {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}


int MSD_limited(float *d_input, float *d_MSD, int nDMs, int nTimesamples, int offset) {
	//---------> Task specific
	ushort nBlocks_x, nBlocks_y;
	int	nBlocks_total, nSteps_x, nSteps_y, nRest;
	float *d_output;
	
	//---------> determining data block size per kernel
	nSteps_x  = PD_NTHREADS;
	nBlocks_x = (int) ((nTimesamples-offset)/nSteps_x);
	nRest     = nTimesamples - offset - nBlocks_x*nSteps_x;
	if(nRest>32) nBlocks_x++;
		
	nSteps_y  = Choose_divider(nDMs,64);
	nBlocks_y = nDMs/nSteps_y;
	nBlocks_total=nBlocks_x*nBlocks_y;
	
	//---------> determining number of threads for final kernel
	int nThreads=2048;
	int itemp=0;
	while(itemp==0 && nThreads>32){
		nThreads=(nThreads>>1);
		itemp=(int) (nBlocks_total/(nThreads*32));
	}
	if(nThreads<32) nThreads=32;
	
	dim3 gridSize(nBlocks_x, nBlocks_y, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(nThreads, 1, 1);
	
	#ifdef MSD_DEBUG
	printf("\n\n");
	printf("----------------> MSD debug: (MSD_limited)\n");
	printf("Kernel for calculating partials:\n");
	printf("ThreadBlocks (TB) in x:%d; Elements processed by TB in x:%d; Remainder in x:%d is processed\n", nBlocks_x, nSteps_x, nRest);
	printf("ThreadBlocks (TB) in y:%d; Elements processed by TB in y:%d; Remainder in y:%d is processed\n", nBlocks_y, nSteps_y, 0);
	printf("gridSize=(%d,%d,%d)\n", gridSize.x, gridSize.y, gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", blockSize.x, blockSize.y, blockSize.z);
	printf("Shared memory required: %0.3f B\n", (float) (PD_NTHREADS*3*4));
	printf("Kernel for final calculation of mean and standard deviation:\n");
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("Memory required for temporary storage:%0.3f MB which is %d floats\n",(nBlocks_total*3*sizeof(float))/(1024.0*1024),  nBlocks_total*3);
	printf("Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	printf("gridSize=(%d,%d,%d)\n", final_gridSize.x, final_gridSize.y, final_gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", final_blockSize.x, final_blockSize.y, final_blockSize.z);	
	printf("---------------------------<\n");
	#endif

	//---------> Allocation of temporary memory
	cudaMalloc((void **) &d_output, nBlocks_total*3*sizeof(float));

	MSD_init();
	MSD_GPU_limited<<<gridSize,blockSize>>>(d_input, d_output, nDMs/nBlocks_y, nTimesamples, offset);
	MSD_GPU_final_regular<<<final_gridSize,final_blockSize>>>(d_output, d_MSD, nBlocks_total);

	cudaFree(d_output);
	
	#ifdef MSD_DEBUG
	float h_MSD[3];
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif

	if (nRest < 32)	return ( nRest );
	else return ( 0 );
}


int MSD_linear_approximation(float *d_input, float *d_MSD_T, int nTaps, int nDMs, int nTimesamples, int offset){
	//---------> Task specific
	int nBlocks_x, nBlocks_y, nBlocks_total, nSteps_x, nSteps_y, nRest, nThreads, itemp; //epw = elements per warp 32 for float 64 for float2
	float *d_output;
	float *d_output_taps;
	float *d_MSD_T_base;

	//---------> determining data block size per kernel
	nSteps_x  = 2*PD_NTHREADS-nTaps+4;
	nBlocks_x = (int) ((nTimesamples-offset)/nSteps_x);
	nRest     = nTimesamples - offset - nBlocks_x*nSteps_x;
	if(nRest>128) nBlocks_x++;
	
	nSteps_y = Choose_divider(nDMs,64);
	nBlocks_y=nDMs/nSteps_y;
	nBlocks_total=nBlocks_x*nBlocks_y;
	
	dim3 gridSize(nBlocks_x, nBlocks_y, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> determining number of threads for final kernel
	nThreads=2048;
	itemp=0;
	while(itemp==0 && nThreads>32){
		nThreads=(nThreads>>1);
		itemp=(int) (nBlocks_total/(nThreads*32));
	}
	if(nThreads<32) nThreads=32;
	
	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(nThreads, 1, 1);

	#ifdef MSD_DEBUG
	printf("\n\n");
	printf("----------------> MSD debug:\n");
	printf("Kernel for calculating partials:\n");
	printf("ThreadBlocks (TB) in x:%d; Elements processed by TB in x:%d; Remainder in x:%d", nBlocks_x, nSteps_x, nRest);
	if(nRest>3*nTaps)//printf(" is processed\n");
	else//printf(" is not processed\n");
	printf("ThreadBlocks (TB) in y:%d; Elements processed by TB in y:%d; Remainder in y:%d is processed\n", nBlocks_y, nSteps_y, 0);
	printf("gridSize=(%d,%d,%d)\n", gridSize.x, gridSize.y, gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", blockSize.x, blockSize.y, blockSize.z);
	printf("Shared memory required: %0.3f B\n", (float) (PD_NTHREADS*24));
	printf("Kernel for final calculation of mean and standard deviation:\n");
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("Memory required for temporary storage:%0.3f MB which is %d floats\n",(nBlocks_total*3*2*sizeof(float))/(1024.0*1024),  nBlocks_total*3*2);
	printf("Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	printf("gridSize=(%d,%d,%d)\n", final_gridSize.x, final_gridSize.y, final_gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", final_blockSize.x, final_blockSize.y, final_blockSize.z);	
	printf("\n");
	#endif
	
	if(nBlocks_total>0){
		//---------> Allocation of temporary memory
		if ( cudaSuccess != cudaMalloc((void **) &d_output, nBlocks_total*3*sizeof(float))) {printf("Allocation error!\n"); exit(1001);}
		if ( cudaSuccess != cudaMalloc((void **) &d_output_taps, nBlocks_total*3*sizeof(float))) {printf("Allocation error!\n"); exit(1001);}
		if ( cudaSuccess != cudaMalloc((void **) &d_MSD_T_base, sizeof(float)*3)) {printf("Allocation error!\n"); exit(1001);}
		
		//---------> MSD
		MSD_init();
		MSD_GPU_LA_ALL<<<gridSize,blockSize>>>(d_input, d_output, d_output_taps, nSteps_y, nTaps, nTimesamples, offset);
		MSD_GPU_final_regular<<<final_gridSize, final_blockSize>>>(d_output, d_MSD_T_base, nBlocks_total);
		MSD_GPU_final_create_LA<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, d_MSD_T_base, nTaps, nBlocks_total);
		
		#ifdef MSD_DEBUG
		float h_MSD_T[3], h_MSD_T_base[3];
		cudaMemcpy(h_MSD_T, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		cudaMemcpy(h_MSD_T_base, d_MSD_T_base, 3*sizeof(float), cudaMemcpyDeviceToHost);
		printf("Output: Mean: %e, Standard deviation: %e; modifier:%e;\n", h_MSD_T[0], h_MSD_T[1], h_MSD_T[2]);
		printf("GPU results after 1 taps: Mean: %e, Standard deviation: %e; Number of elements:%d;\n", h_MSD_T_base[0], h_MSD_T_base[1], (int) h_MSD_T_base[2]);
		printf("---------------------------<\n");
		#endif
		
		//---------> De-allocation of temporary memory
		cudaFree(d_output);
		cudaFree(d_output_taps);
		cudaFree(d_MSD_T_base);
	}
	else {
		printf("Number of time samples is too small! Increase number of samples send to the boxcar filters. (MSD_linear_approximation)\n");
		exit(1002);
	}
	
	if(nRest<64) return(nRest);
	else return(0);	
}


int MSD_LA_Nth(float *d_input, float *d_bv_in, float *d_MSD_T, float *d_MSD_DIT, int nTaps, int nDMs, int nTimesamples, int offset, int DIT_value){
	//---------> Task specific
	int nBlocks_x, nBlocks_y, nBlocks_total, nSteps_x, nSteps_y, nRest, nThreads, itemp; //epw = elements per warp 32 for float 64 for float2
	float *d_output;
	float *d_output_FIR;
	float *d_MSD_T_base;

	//---------> determining data block size per kernel
	nSteps_x  = 2*PD_NTHREADS-nTaps+4;
	nBlocks_x = (int) ((nTimesamples-offset)/nSteps_x);
	nRest     = nTimesamples - offset - nBlocks_x*nSteps_x;
	if(nRest>0) nBlocks_x++;
	
	nSteps_y = Choose_divider(nDMs,64);
	nBlocks_y=nDMs/nSteps_y;
	nBlocks_total=nBlocks_x*nBlocks_y;
	
	dim3 gridSize(nBlocks_x, nBlocks_y, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> determining number of threads for final kernel
	nThreads=2048;
	itemp=0;
	while(itemp==0 && nThreads>32){
		nThreads=(nThreads>>1);
		itemp=(int) (nBlocks_total/(nThreads*32));
	}
	if(nThreads<32) nThreads=32;
	
	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(nThreads, 1, 1);
	
	#ifdef MSD_DEBUG
	printf("\n\n");
	printf("----------------> MSD debug:\n");
	printf("Kernel for calculating partials: (MSD_LA_Nth)\n");
	printf("nTimesamples:%d; offset:%d, nDMs:%d;\n", nTimesamples, offset, nDMs);
	printf("ThreadBlocks (TB) in x:%d; Elements processed by TB in x:%d; Remainder in x:%d", nBlocks_x, nSteps_x, nRest);
	if(nRest>3*nTaps)//printf(" is processed\n");
	else//printf(" is not processed\n");
	printf("ThreadBlocks (TB) in y:%d; Elements processed by TB in y:%d; Remainder in y:%d is processed\n", nBlocks_y, nSteps_y, 0);
	printf("gridSize=(%d,%d,%d)\n", gridSize.x, gridSize.y, gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", blockSize.x, blockSize.y, blockSize.z);
	printf("Shared memory required: %0.3f B\n", (float) (PD_NTHREADS*24));
	printf("Kernel for final calculation of mean and standard deviation:\n");
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("Memory required for temporary storage:%0.3f MB which is %d floats\n",(nBlocks_total*3*2*sizeof(float))/(1024.0*1024),  nBlocks_total*3*2);
	printf("Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	printf("gridSize=(%d,%d,%d)\n", final_gridSize.x, final_gridSize.y, final_gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", final_blockSize.x, final_blockSize.y, final_blockSize.z);	
	printf("\n");
	#endif
	
	if(nBlocks_total>0){
		//---------> Allocation of temporary memory
		if ( cudaSuccess != cudaMalloc((void **) &d_output, nBlocks_total*3*sizeof(float))) {printf("Allocation error!\n"); exit(1001);}
		if ( cudaSuccess != cudaMalloc((void **) &d_output_FIR, nBlocks_total*3*sizeof(float))) {printf("Allocation error!\n"); exit(1001);}
		if ( cudaSuccess != cudaMalloc((void **) &d_MSD_T_base, sizeof(float)*3)) {printf("Allocation error!\n"); exit(1001);}
		
		//---------> MSD
		MSD_init();
		MSD_GPU_LA_ALL_Nth<<<gridSize,blockSize>>>(d_input, d_bv_in, d_output, d_output_FIR, nSteps_y, nTaps, nTimesamples, offset);
		MSD_GPU_final_regular<<<final_gridSize, final_blockSize>>>(d_output, d_MSD_T_base, nBlocks_total);
		MSD_GPU_final_create_LA_Nth<<<final_gridSize, final_blockSize>>>(d_output_FIR, d_MSD_T, d_MSD_T_base, d_MSD_DIT, nTaps, nBlocks_total, DIT_value);
		
		#ifdef MSD_DEBUG
		float h_MSD_T[4], h_MSD_T_base[3], h_MSD_DIT[3];
		cudaMemcpy(h_MSD_T, d_MSD_T, 4*sizeof(float), cudaMemcpyDeviceToHost); 
		cudaMemcpy(h_MSD_T_base, d_MSD_T_base, 3*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(h_MSD_DIT, d_MSD_DIT, 3*sizeof(float), cudaMemcpyDeviceToHost);
		printf("d_MSD_T: BV Mean: %f, Standard deviation: %f; modifier:%f; DIT Mean:%f\n", h_MSD_T[0], h_MSD_T[1], h_MSD_T[2], h_MSD_T[3]);
		printf("MSD for d_bv_in: Mean: %f, Standard deviation: %f; Number of elements:%d;\n", h_MSD_T_base[0], h_MSD_T_base[1], (int) h_MSD_T_base[2]);
		printf("MSD for DIT: Mean: %f, Standard deviation: %f; Number of elements:%d;\n", h_MSD_DIT[0], h_MSD_DIT[1], (int) h_MSD_DIT[2]);
		printf("---------------------------<\n");
		#endif
		
		//---------> De-allocation of temporary memory
		cudaFree(d_output);
		cudaFree(d_output_FIR);
		cudaFree(d_MSD_T_base);
	}
	else {
		printf("WARNING: Number of time samples is too small! Increase number of samples send to the boxcar filters. (MSD_LA_Nth)\n");
		return(1);
	}

	return(0);
}

