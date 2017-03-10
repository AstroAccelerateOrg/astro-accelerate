//Added by Karel Adamek

#include "headers/params.h"
#include "device_MSD_plane_kernel.cu"

void MSD_init(void)
{
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}

int MSD(float *d_input, float *d_MSD, int nDMs, int nTimesamples)
{
	//---------> Task specific
	int nBlocks, nRest, Elements_per_block, nElements, nThreads;
	float *d_output;

	nElements = nDMs * nTimesamples;
	Elements_per_block = 2 * WARP * MSD_ELEM_PER_THREAD * MSD_WARPS_PER_BLOCK;
	nBlocks = nElements / Elements_per_block;
	nRest = nElements - nBlocks * Elements_per_block;

	//---------> CUDA block and CUDA grid parameters
	int nCUDAblocks_x = nBlocks;
	int nCUDAblocks_y = 1;

	dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);
	dim3 blockSize(WARP * MSD_WARPS_PER_BLOCK, 1, 1);

	if (nRest < 128)
		nThreads = 32;
	else
		nThreads = 128;

	dim3 remainder_gridSize(1, 1, 1);
	dim3 remainder_blockSize(nThreads, 1, 1);

	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(WARP * MSD_WARPS_PER_BLOCK, 1, 1);

	//---------> Allocation of temporary memory
	cudaMalloc((void **) &d_output, ( nBlocks * 2 + 3 ) * sizeof(float));

	//---------> Pulse detection FIR
	MSD_init();
	MSD_GPU<<<gridSize, blockSize>>>((float2 *) d_input, d_output);
	if (nRest > 0)
	{
		MSD_GPU_remainder<<<remainder_gridSize, remainder_blockSize>>>(&d_input[nBlocks * 2 * WARP * MSD_ELEM_PER_THREAD * MSD_WARPS_PER_BLOCK], &d_output[2 * nBlocks], nRest);
	}
	MSD_GPU_final<<<final_gridSize, final_blockSize>>>(d_output, d_MSD, nBlocks, nRest, (float) nElements);

	//---------> De-allocation of temporary memory
	cudaFree(d_output);

	// Unprocessed samples depends on whether kernel MSD_GPU_one_reduction_remainder is launched or not. If not then unprocessed samples are nRest BUT only at the end since this treats data as 1d array! 
	return ( 0 );
}

