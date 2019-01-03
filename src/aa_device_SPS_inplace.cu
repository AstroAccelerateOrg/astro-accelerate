//Added by Karel Adamek

#include "aa_device_SPS_inplace.hpp"

namespace astroaccelerate {

  int Choose_dim_SPS(int grid_dim)
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

  void PD_SEARCH_INPLACE_init()
  {
    //---------> Specific nVidia stuff
    cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
    cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeEightByte);
  }

  int PD_SEARCH_INPLACE(float *d_input, unsigned char *d_output_taps, float *d_MSD, int maxTaps, int nDMs, int nTimesamples)
  {
    //---------> Task specific
    int nBlocks_x, nBlocks_y, nRest, Elements_per_block, nLoops, nThreads;
    int SM_size = ( PD_NTHREADS * PD_NWINDOWS + maxTaps - 1 ) * 4;
    float *d_output;

    Elements_per_block = PD_NTHREADS * PD_NWINDOWS;
    nBlocks_x = nTimesamples / Elements_per_block;
    nRest = nTimesamples - nBlocks_x * Elements_per_block;
    if (nRest > 0)
      nBlocks_x++;

    nLoops = maxTaps / WARP;
    nRest = maxTaps - nLoops * WARP;
    if (nRest > 0)
      nLoops++;

    nThreads = Choose_dim_SPS(nDMs);
    nBlocks_y = nDMs / nThreads;

    cudaMalloc((void **) &d_output, nBlocks_x * ( maxTaps - 1 ) * nDMs * sizeof(float));
    //cudaMemset((void*) d_output, 0, nBlocks_x*(maxTaps-1)*nDMs*sizeof(float));

    PD_SEARCH_INPLACE_init();

    //---------> CUDA block and CUDA grid parameters for temporary storage
    dim3 gridSize_temp(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize_temp(WARP, nThreads, 1);
    call_kernel_PD_ZC_GPU_KERNEL(gridSize_temp, blockSize_temp, d_input, d_output, maxTaps,
				 nTimesamples, nLoops);

    //---------> CUDA block and CUDA grid parameters for in-place PD
    dim3 gridSize(nBlocks_x, nDMs, 1);
    dim3 blockSize(PD_NTHREADS, 1, 1);
    call_kernel_PD_INPLACE_GPU_KERNEL(gridSize, blockSize, SM_size, d_input, d_output, d_output_taps,
				      d_MSD, maxTaps, nTimesamples);

    cudaFree(d_output);

    return ( maxTaps );
  }

} //namespace astroaccelerate
