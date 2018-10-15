#ifndef ASTRO_ACCELERATE_DEVICE_POWER_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_POWER_KERNEL_HPP

#include <cufft.h>

void call_kernel_power_kernel(const dim3 &        block_size,
                              const dim3 &        grid_size,
                              const int &         smem_bytes,
                              const cudaStream_t &stream,
                              const int &         half_samps,
                              const int &         acc,
                              cufftComplex *const d_signal_fft,
                              float *const        d_signal_power);
void call_kernel_GPU_simple_power_and_interbin_kernel(
    const dim3 &  grid_size,
    const dim3 &  block_size,
    float2 *const d_input_complex,
    float *const  d_output_power,
    float *const  d_output_interbinning,
    const int &   nTimesamples,
    const float & norm);

#endif
