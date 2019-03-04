#include <stdio.h>
#include "aa_params.hpp"
#include "aa_device_set_stretch.hpp"
#include "aa_device_set_stretch_kernel.hpp"

namespace astroaccelerate {
  /** \brief Function for Doppler stretch. */
  void set_stretch_gpu(cudaEvent_t event, cudaStream_t stream, int samps, float mean, float *d_input) {

    int divisions_in_t = 32;
    int num_blocks_t = samps / divisions_in_t;

    dim3 threads_per_block(divisions_in_t);
    dim3 num_blocks(num_blocks_t);

    cudaStreamWaitEvent(stream, event, 0);
    call_kernel_set_stretch_kernel(num_blocks, threads_per_block, 0, stream, samps, mean, d_input);
    //getLastCudaError("stretch_kernel failed");
    cudaEventRecord(event, stream);
  }
} //namespace astroaccelerate
