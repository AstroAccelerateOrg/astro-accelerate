#include "aa_device_SNR_limited.hpp"

namespace astroaccelerate {

  int Choose_dim(int grid_dim) {
    int seive[15] =
      { 32, 31, 29, 23, 19, 17, 16, 13, 11, 8, 7, 5, 4, 3, 2 };

    int f, nRest, nBlocks, N, N_accepted;

    N = 1;
    N_accepted = 1;
    for (int i = 0; i < 4; i++) {
      for (f = 0; f < 15; f++) {
	nBlocks = grid_dim / seive[f];
	nRest = grid_dim - nBlocks * seive[f];
	if (nRest == 0) {
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

  void SNR_limited_init() {
    //---------> Specific nVidia stuff
    cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
    cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
  }

  int SNR_limited(float *d_FIR_input, float *d_SNR_output, ushort *d_SNR_taps, float *d_MSD, int nTaps, int nDMs, int nTimesamples, int offset) {
    //---------> Task specific
    int nBlocks_x, nBlocks_y, nSteps_x, nSteps_y, nRest, nThreads, // int nBlocks_total;
      epw; //epw = elements per warp 32 for float 64 for float2
    // int nElements;
    //---------> CUDA block and CUDA grid parameters
    // Determining in x direction (direction of data alignment)
    epw = 32;
    nBlocks_x = 0;
    nRest = 0;

    nSteps_x = Choose_dim(( nTimesamples ) / epw);
    nBlocks_x = nBlocks_x + ( nTimesamples - offset ) / ( nSteps_x * epw );
    nRest += nTimesamples - offset - nBlocks_x * nSteps_x * epw;
    if (nRest > 0)
      nBlocks_x++;

    nSteps_y = Choose_dim(nDMs);
    nBlocks_y = nDMs / nSteps_y;

    //nBlocks_total = nBlocks_x * nBlocks_y;
    //nElements = nBlocks_total * nSteps_x * epw * nSteps_y;

    nThreads = nSteps_y * WARP;

    // calculation of the partials
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(nThreads, 1, 1);

    //---------> MSD
    SNR_limited_init();
    call_kernel_SNR_GPU_limited(gridSize, blockSize, d_FIR_input, d_SNR_output, d_SNR_taps, d_MSD,
				nSteps_x, nTaps, nTimesamples, offset);

    if (nRest < epw)
      return ( nRest );
    else
      return ( 0 );

  }

} //namespace astroaccelerate
