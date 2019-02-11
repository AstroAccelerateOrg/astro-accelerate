//Added by Karel Adamek

#include <stdio.h>
#include "aa_params.hpp"
#include "aa_device_single_pulse_search_kernel.hpp"

namespace astroaccelerate {

  void PD_SEARCH_init(void)
  {
    //---------> Specific nVidia stuff
    cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
    cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeEightByte);
  }

  int PD_SEARCH(float *d_input, float *d_output, float *d_output_taps, float *d_MSD, int maxTaps, int nDMs, int nTimesamples)
  {
    //---------> Task specific
    int nBlocks, nRest, Elements_per_block, ut;

    //--------->  Moved to device_load_data.cu
    //float h_sqrt_taps[PD_MAXTAPS+1];
    //for(int f=0; f<=PD_MAXTAPS; f++) h_sqrt_taps[f]=(float) sqrt((double) f);
    //cudaMemcpyToSymbol(c_sqrt_taps, h_sqrt_taps, (PD_MAXTAPS+1)*sizeof(float));

    Elements_per_block = PD_NTHREADS * PD_NWINDOWS;
    nBlocks = ( nTimesamples - maxTaps + 1 ) / Elements_per_block;
    nRest = ( nTimesamples - maxTaps + 1 ) - nBlocks * Elements_per_block;
    if (nRest > 0)
      printf("nRest:%d\n", nRest);

    //---------> CUDA block and CUDA grid parameters
    int nCUDAblocks_x = nBlocks;
    int nCUDAblocks_y = nDMs;
    int SM_size = ( PD_NTHREADS * PD_NWINDOWS + maxTaps - 1 ) * 4;

    dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);
    dim3 blockSize(PD_NTHREADS, 1, 1);

    //---------> Pulse detection FIR
    PD_SEARCH_init();
    call_kernel_PD_SEARCH_GPU(gridSize, blockSize, SM_size, d_input, d_output, d_output_taps, d_MSD, maxTaps, nTimesamples);

    ut = nRest + maxTaps - 1;
    return ( ut );
  }

} //namespace astroaccelerate
