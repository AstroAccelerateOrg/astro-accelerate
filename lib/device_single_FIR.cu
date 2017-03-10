//Added by Karel Adamek

#include "headers/params.h"

void PD_FIR_init(void)
{
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeEightByte);
}

int PD_FIR(float *d_input, float *d_output, int nTaps, int nDMs, int nTimesamples)
{
	//---------> Task specific
	int ut; //unused timesamples
	int itemp = (int) ( ( nTaps - 1 ) / ( WARP * PD_FIR_ACTIVE_WARPS ) ) + 1;
	int nLoops = PD_FIR_NWINDOWS + itemp;

	//---------> CUDA block and CUDA grid parameters
	int nCUDAblocks_x = (int) ( ( nTimesamples - nTaps + 1 ) / ( PD_FIR_ACTIVE_WARPS * WARP * PD_FIR_NWINDOWS ) );
	int nCUDAblocks_y = nDMs;
	int SM_size = ( PD_FIR_ACTIVE_WARPS * WARP * PD_FIR_NWINDOWS + nTaps - 1 ) * 4;

	dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);
	dim3 blockSize(PD_FIR_ACTIVE_WARPS * WARP, 1, 1);

	//---------> Pulse detection FIR
	PD_FIR_init();
	PD_FIR_GPU<<<gridSize, blockSize, SM_size>>>(d_input, d_output, nTaps, nLoops, nTimesamples);

	ut = nTimesamples - nCUDAblocks_x * PD_FIR_ACTIVE_WARPS * WARP * PD_FIR_NWINDOWS;
	return ( ut );
}
