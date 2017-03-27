//Added by Karel Adamek

#include "headers/params.h"
#include "device_MSD_limited_kernel.cu"

int Choose_x_dim(int grid_dim)
{
	int seive[15] =
	{ 32, 31, 29, 23, 19, 17, 16, 13, 11, 8, 7, 5, 4, 3, 2 };

	int f, nRest, nBlocks, N, N_accepted;

	N = 1;
	N_accepted = 1;
	for (int i = 0; i < 4; i++)
	{
		for (f = 0; f < 15; f++)
		{
			nBlocks = grid_dim / seive[f];
			nRest = grid_dim - nBlocks * seive[f];
			if (nRest == 0)
			{
				N_accepted = N_accepted * N;
				N = seive[f];
				break;
			}
		}
		if (( N_accepted * N ) > 32 || N == 1)
			return ( N_accepted );
		grid_dim = grid_dim / N;
	}

	return ( N_accepted );
}

int Choose_y_dim(int grid_dim)
{
	int seive[5] =
	{ 32, 16, 8, 4, 2 };

	int f, nRest, nBlocks, N;

	N = 1;
	for (f = 0; f < 5; f++)
	{
		nBlocks = grid_dim / seive[f];
		nRest = grid_dim - nBlocks * seive[f];
		if (nRest == 0)
		{
			N = seive[f];
			break;
		}
	}

	return ( N );
}

void MSD_limited_init()
{
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}

int MSD_limited(float *d_input, float *d_MSD, int nDMs, int nTimesamples, int offset)
{
	//---------> Task specific
	int nBlocks_x, nBlocks_y, nBlocks_total, nSteps_x, nSteps_y, nRest, nThreads,
	    epw; //epw = elements per warp 32 for float 64 for float2
	// int nElements;
	float *d_output;

	//---------> CUDA block and CUDA grid parameters
	// Determining in x direction (direction of data alignment)
	epw = 32;
	nBlocks_x = 0;
	nRest = 0;
	// nTimesamples must be divisible by 64 because otherwise it is extremely tedious and I refuse to do it.
	nSteps_x = Choose_x_dim(( nTimesamples ) / epw);
	nBlocks_x = nBlocks_x + ( nTimesamples - offset ) / ( nSteps_x * epw );
	nRest += nTimesamples - offset - nBlocks_x * nSteps_x * epw;
	if (nRest > epw)
		nBlocks_x++; // if nRest<64 then it means a lot of branching in the kernel and error it induces would be generally small.

	nSteps_y = Choose_y_dim(nDMs);
	nBlocks_y = nDMs / nSteps_y;
	// I do not calculate nRest here since I assume it will be always divisible by nSteps_y.

	nBlocks_total = nBlocks_x * nBlocks_y;
	//nElements = nBlocks_total * nSteps_x * epw * nSteps_y;

	nThreads = nSteps_y * WARP;

	// calculation of the partials
	dim3 gridSize(nBlocks_x, nBlocks_y, 1);
	dim3 blockSize(nThreads, 1, 1);

	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(WARP * 4, 1, 1);

	//---------> Allocation of temporary memory
	cudaMalloc((void **) &d_output, nBlocks_total * 3 * sizeof(float));

	//---------> MSD
	MSD_init();
	MSD_GPU_limited<<<gridSize, blockSize, nThreads * 12>>>(d_input, d_output, nSteps_x, nTimesamples, offset);
	MSD_GPU_limited_final<<<final_gridSize, final_blockSize>>>(d_output, d_MSD, nBlocks_total);

	//---------> De-allocation of temporary memory
	cudaFree(d_output);

	if (nRest < 64)
		return ( nRest );
	else
		return ( 0 );
}

