#include <stdio.h>
#include "aa_params.hpp"
#include "aa_device_stretch_kernel.hpp"

namespace astroaccelerate {
  /** \brief Doppler stretch. Simple corner turn on the GPU. */
  void stretch_gpu(cudaEvent_t event, cudaStream_t stream, int acc, int samps, float tsamp, float *d_input, float *d_output) {
    int divisions_in_t = 32;
    int num_blocks_t = samps / divisions_in_t;

    float t_zero = ( (double) tsamp ) / ( 1.0 + ( ( acc * samps * (double) tsamp ) / 599584916.0 ) );
    float multiplier = ( t_zero * acc * (double) tsamp ) / 599584916.0;
    float tsamp_inverse = 1.0 / tsamp;

    dim3 threads_per_block(divisions_in_t);
    dim3 num_blocks(num_blocks_t);

    cudaStreamWaitEvent(stream, event, 0);
    call_kernel_stretch_kernel(num_blocks, threads_per_block, 0, stream, acc, samps, tsamp, d_input, d_output, t_zero, multiplier, tsamp_inverse);
    //getLastCudaError("stretch_kernel failed");
    cudaEventRecord(event, stream);
  }
} //namespace astroaccelerate
