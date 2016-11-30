//Added by Karel Adamek

#include "AstroAccelerate/params.h"
#include "device_threshold_kernel.cu"

void THR_init(void)
{
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}

int THRESHOLD(float *d_input, unsigned char *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nDMs, int nTimesamples, int offset, int max_list_size)
{
	//---------> Task specific
	int nBlocks, nRest, Elements_per_block;

	Elements_per_block = 2*WARP*THR_ELEM_PER_THREAD;
	nBlocks = (nTimesamples-offset)/Elements_per_block;
	nRest = (nTimesamples-offset) - nBlocks*Elements_per_block;
	if(nRest>0) nBlocks++;

	//---------> CUDA block and CUDA grid parameters
	int nCUDAblocks_x = nBlocks;
	int nCUDAblocks_y = nDMs/THR_WARPS_PER_BLOCK;

	dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);
	dim3 blockSize(WARP*THR_WARPS_PER_BLOCK, 1, 1);

	//---------> Pulse detection FIR
	THR_init();
	THR_GPU_WARP<<<gridSize, blockSize>>>((float2 *) d_input, d_input_taps, d_output_list, gmem_pos, threshold, nTimesamples, nTimesamples-offset, max_list_size);

	return ( 0 );
}

int THRESHOLD_ignore(float *d_input, unsigned char *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nTaps, int nDMs, int nTimesamples, int offset, int max_list_size)
{
	//---------> Task specific
	int nBlocks, nRest, Elements_per_block;

	Elements_per_block = 2*WARP*THR_ELEM_PER_THREAD;
	nBlocks = (nTimesamples-offset)/Elements_per_block;
	nRest = nTimesamples - nBlocks*Elements_per_block;

	//---------> CUDA block and CUDA grid parameters
	int nCUDAblocks_x = nBlocks;
	int nCUDAblocks_y = nDMs/THR_WARPS_PER_BLOCK;

	dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);
	dim3 blockSize(WARP*THR_WARPS_PER_BLOCK, 1, 1);

	//---------> Pulse detection FIR
	THR_init();
	THR_GPU_WARP_ignore<<<gridSize, blockSize>>>((float2 *) d_input, d_input_taps, d_output_list, gmem_pos, threshold, nTaps, nTimesamples, max_list_size);

	return ( 0 );
}
