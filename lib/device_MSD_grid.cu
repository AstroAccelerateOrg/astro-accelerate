//Added by Karel Adamek
#include "headers/params.h"
#include "device_MSD_grid_kernel.cu"

void MSD_grid_init(void)
{
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferNone);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}

int MSD_grid(float *d_input, float *d_output, int CellDim_x, int CellDim_y, int nDMs, int nTimesamples)
{
	//---------> Task specific
	int GridSize_x, GridSize_y, x_steps, y_steps, nThreads;
	GridSize_x = nTimesamples / CellDim_x;
	GridSize_y = nDMs / CellDim_y;
	x_steps = CellDim_x / WARP;
	if (CellDim_y < 16)
	{
		y_steps = 1;
		nThreads = WARP * CellDim_y;
	}
	else
	{
		nThreads = WARP * 16;
		y_steps = CellDim_y / 16;
	}
	//---------> Initial phase
	int nCUDAblocks_x = GridSize_x;
	int nCUDAblocks_y = GridSize_y;

	dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);
	dim3 blockSize(nThreads, 1, 1);

	//---------> Pulse detection FIR
	MSD_grid_init();
	MSD_GPU_grid<<<gridSize, blockSize, nThreads * 8>>>(d_input, d_output, x_steps, y_steps, nTimesamples, 1);

	return ( 1 );
}
