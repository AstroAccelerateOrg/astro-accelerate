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

void MSD_limited_init() {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}


int MSD_limited(float *d_input, float *d_MSD, int nDMs, int nTimesamples, int offset) {
	//---------> Task specific
	int nBlocks_x, nBlocks_y, nBlocks_total, nSteps_x, nSteps_y, nRest, nThreads, epw; //epw = elements per warp 32 for float 64 for float2
	float *d_output;

	//---------> CUDA block and CUDA grid parameters
	// Determining in x direction (direction of data alignment)
	epw = 32;
	nBlocks_x = 0;
	nRest = 0;
	
	nSteps_x = Choose_x_dim(( nTimesamples ) / epw);
	nBlocks_x = nBlocks_x + ( nTimesamples - offset ) / ( nSteps_x*epw );
	nRest += nTimesamples - offset - nBlocks_x*nSteps_x*epw;
	if (nRest > epw) nBlocks_x++; // if nRest<64(32) then it means a lot of branching in the kernel and error it induces would be generally small.

	nSteps_y = Choose_y_dim(nDMs);
	nBlocks_y = nDMs / nSteps_y;
	// I do not calculate nRest here since I assume it will be always divisible by nSteps_y.

	nBlocks_total = nBlocks_x*nBlocks_y;
	//nElements = nBlocks_total*nSteps_x*epw*nSteps_y;

	nThreads = nSteps_y*WARP;

	// calculation of the partials
	dim3 gridSize(nBlocks_x, nBlocks_y, 1);
	dim3 blockSize(nThreads, 1, 1);

	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(WARP*4, 1, 1);
	
	#ifdef MSD_DEBUG
	printf("\n\n");
	printf("----------------> MSD debug:\n");
	printf("Kernel for calculating partials:\n");
	printf("ThreadBlocks (TB) in x:%d; Elements processed by TB in x:%d; Remainder in x:%d is processed\n", nBlocks_x, nSteps_x, nRest);
	printf("ThreadBlocks (TB) in y:%d; Elements processed by TB in y:%d; Remainder in y:%d is processed\n", nBlocks_y, nSteps_y, 0);
	printf("gridSize=(%d,%d,%d)\n", gridSize.x, gridSize.y, gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", blockSize.x, blockSize.y, blockSize.z);
	printf("Shared memory required: %0.3f B\n", (float) (nThreads*12));
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

	//---------> MSD
	MSD_init();
	MSD_GPU_limited<<<gridSize, blockSize, nThreads*12>>>(d_input, d_output, nSteps_x, nTimesamples, offset);
	MSD_GPU_limited_final<<<final_gridSize, final_blockSize>>>(d_output, d_MSD, nBlocks_total);

	//---------> De-allocation of temporary memory
	cudaFree(d_output);

	if (nRest < 64)
		return ( nRest );
	else
		return ( 0 );
}

int MSD_linear_approximation(float *d_input, float *d_MSD_T, int nTaps, int nDMs, int nTimesamples, int offset){
	//---------> Task specific
	int nBlocks_x, nBlocks_y, nBlocks_total, nSteps_x, nSteps_y, nRest, nThreads, itemp; //epw = elements per warp 32 for float 64 for float2
	float *d_output;
	float *d_output_taps;
	float *d_MSD_T_base;

	//---------> CUDA block and CUDA grid parameters
	// Determining in x direction (direction of data alignment)
	
	nBlocks_x=0; nRest=0;

	nSteps_x=2*PD_NTHREADS-nTaps;
	nBlocks_x=nBlocks_x + (nTimesamples-offset)/(nSteps_x);
	nRest=nTimesamples - offset - nBlocks_x*nSteps_x;
	if(nRest>0) nBlocks_x++; // if nRest<64 then it means a lot of branching in the kernel and error it induces would be generally small.
	//printf("nSteps_x:%d;b nBlocks_x:%d; nRest:%d; \n", nSteps_x, nBlocks_x, nRest);
	
	
	nSteps_y=Choose_y_dim(nDMs);
	nBlocks_y=nDMs/nSteps_y;
	// I do not calculate nRest here since I assume it will be always divisible by nSteps_y.
	//printf("nSteps_y:%d; nBlocks_y:%d; nRest:%d; \n", nSteps_y, nBlocks_y, 0);
	
	nBlocks_total=nBlocks_x*nBlocks_y;
	
	// calculation of the partials
	dim3 gridSize(nBlocks_x, nBlocks_y, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
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
	printf("ThreadBlocks (TB) in x:%d; Elements processed by TB in x:%d; Remainder in x:%d is processed\n", nBlocks_x, nSteps_x, nRest);
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
	printf("---------------------------<\n");
	#endif
	
	// ----------------------------------------------->
	// --------> Measured part (Pulse detection FIR)
	
	//---------> Allocation of temporary memory
	if ( cudaSuccess != cudaMalloc((void **) &d_output, nBlocks_total*3*sizeof(float))) {printf("Allocation error!\n"); exit(1001);}
	if ( cudaSuccess != cudaMalloc((void **) &d_output_taps, nBlocks_total*3*sizeof(float))) {printf("Allocation error!\n"); exit(1001);}
	if ( cudaSuccess != cudaMalloc((void **) &d_MSD_T_base, sizeof(float)*3)) {printf("Allocation error!\n"); exit(1001);}
	
	
	//---------> MSD
	MSD_init();
	MSD_GPU_LA_ALL<<<gridSize,blockSize>>>(d_input, d_output, d_output_taps, nSteps_y, nTaps, nTimesamples, offset);
	MSD_GPU_limited_final<<<final_gridSize, final_blockSize>>>(d_output, d_MSD_T_base, nBlocks_total);
	MSD_GPU_limited_final_create_linear_approx<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, d_MSD_T_base, nTaps, nBlocks_total);
		
	//---------> De-allocation of temporary memory
	cudaFree(d_output);
	cudaFree(d_output_taps);
	cudaFree(d_MSD_T_base);
	
	// --------> Measured part (Pulse detection FIR)
	// ----------------------------------------------->
	
	if(nRest<64) return(nRest);
	else return(0);
	
}

