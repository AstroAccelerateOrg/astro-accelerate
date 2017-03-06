//Added by Karel Adamek

#include "AstroAccelerate/params.h"
#include "device_BLN_kernel.cu"

void BLN_init(){
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}

int BLN(float *d_input, float *d_MSD, int CellDim_x, int CellDim_y, int nDMs, int nTimesamples, int offset, float multiplier){
	//---------> Task specific
	int GridSize_x, GridSize_y, x_steps, y_steps, nThreads;
	GridSize_x=(nTimesamples-offset)/CellDim_x;
	GridSize_y=nDMs/CellDim_y;
	x_steps=CellDim_x/WARP;
	if(CellDim_y<HALF_WARP) {
		y_steps  = 1;
		nThreads = WARP*CellDim_y;
	}
	else {
		nThreads = WARP*HALF_WARP;
		y_steps  = CellDim_y/HALF_WARP;
	}

	//---------> Initial phase
	dim3 gridSize(GridSize_x, GridSize_y, 1);
	dim3 blockSize(nThreads, 1, 1);

	//---------> Final phase
	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(WARP*WARP, 1, 1);

	//---------> Allocation of temporary memory
	float *d_output;
	cudaMalloc((void **) &d_output, GridSize_x*GridSize_y*3*sizeof(float));

	//---------> MSD
	BLN_init();
	BLN_MSD_GPU_grid<<<gridSize,blockSize,nThreads*8>>>(d_input, d_output, x_steps, y_steps, nTimesamples, 0);
	BLN_outlier_rejection<<<final_gridSize, final_blockSize>>>(d_output, d_MSD, GridSize_x*GridSize_y, (float) (CellDim_x*CellDim_y), multiplier);

	//---------> De-allocation of temporary memory
	cudaFree(d_output);
	
	return(1);
}

