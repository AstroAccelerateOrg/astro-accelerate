//Added by Karel Adamek
//#define MSD_DEBUG

#include "headers/params.h"
#include "headers/device_MSD_Configuration.h"
#include "device_MSD_limited_kernel.cu"

void MSD_init(void) {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
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

/*
void MSD_limited_init() {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}
*/

/*
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
*/

int MSD_limited(float *d_input, float *d_MSD, float *d_temp, MSD_Configuration *MSD_conf) {
	
	#ifdef MSD_DEBUG
	MSD_conf->print();
	#endif

	MSD_init();
	MSD_GPU_limited<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
	
	#ifdef MSD_DEBUG
	float h_MSD[3];
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif

	return (0);
}

int MSD_limited(float *d_input, float *d_MSD, int nDMs, int nTimesamples, int offset){
	int result;
	MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
	float *d_temp;
	cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));
	result = MSD_limited(d_input, d_MSD, d_temp, &conf);
	cudaFree(d_temp);
	return(result);
}


/*
int MSD_limited_continuous(float *d_input, float *d_MSD, float *d_previous_partials, int nDMs, int nTimesamples, int offset) {
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
	MSD_GPU_final_regular<<<final_gridSize,final_blockSize>>>(d_output, d_MSD, d_previous_partials, nBlocks_total);

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
*/



int MSD_limited_continuous(float *d_input, float *d_MSD, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf) {

	#ifdef MSD_DEBUG
	MSD_conf->print();
	#endif

	MSD_init();
	MSD_GPU_limited<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, d_previous_partials, MSD_conf->nBlocks_total);
	
	#ifdef MSD_DEBUG
	float h_MSD[3];
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif

	return (0);
}

int MSD_limited_continuous(float *d_input, float *d_MSD, float *d_previous_partials, int nDMs, int nTimesamples, int offset){
	int result;
	MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
	float *d_temp;
	cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));
	result = MSD_limited_continuous(d_input, d_MSD, d_previous_partials, d_temp, &conf);
	cudaFree(d_temp);
	return(result);
}

