#include <stdio.h>
#include <cufft.h>
#include "aa_params.hpp"
#include "aa_device_power_kernel.hpp"

// define for debug info
//#define POWER_DEBUG

//{{{ Dopler Stretch 

namespace astroaccelerate {

  void power_gpu(cudaEvent_t event, cudaStream_t stream, int samps, int acc, cufftComplex *d_signal_fft, float *d_signal_power)
  {

    int half_samps = samps / 2;

    int divisions_in_t = 32;
    int num_blocks_t = half_samps / divisions_in_t;

    dim3 threads_per_block(divisions_in_t);
    dim3 num_blocks(num_blocks_t);

    cudaStreamWaitEvent(stream, event, 0);
    call_kernel_power_kernel(num_blocks, threads_per_block, 0, stream, half_samps, acc, d_signal_fft, d_signal_power);
    //getLastCudaError("power_kernel failed");
    cudaEventRecord(event, stream);
  }

  //}}}


  void simple_power_and_interbin(float2 *d_input, float *d_power_output, float *d_interbin_output, int nTimesamples, int nDMs){
    int n_blocks_x, n_blocks_y, itemp;
	
    n_blocks_x = (nTimesamples>>1)/PAI_NTHREADS;
    itemp = (nTimesamples>>1) - n_blocks_x*PAI_NTHREADS;
    if (itemp>0) n_blocks_x++;
    n_blocks_y = nDMs;
	
    dim3 blockDim(PAI_NTHREADS, 1, 1);
    dim3 gridSize(n_blocks_x ,n_blocks_y , 1);
	
#ifdef POWER_DEBUG
    printf("gridSize=(%d,%d,%d)\n", gridSize.x, gridSize.y, gridSize.z);
    printf("blockDim=(%d,%d,%d)\n", blockDim.x, blockDim.y, blockDim.z);
#endif
	
    call_kernel_GPU_simple_power_and_interbin_kernel(gridSize,blockDim, d_input, d_power_output, d_interbin_output, nTimesamples, sqrt(nTimesamples));
  }

} //namespace astroaccelerate
