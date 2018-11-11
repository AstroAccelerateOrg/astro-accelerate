#ifndef ASTRO_ACCELERATE_DEVICE_SNR_LIMITED_INPLACE_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_SNR_LIMITED_INPLACE_KERNEL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include "params.hpp"
#include "device_SPS_inplace_kernel.hpp"

namespace astroaccelerate {

void call_kernel_PD_ZC_GPU(float *const d_input, float *const d_output, const int &maxTaps, const int &nTimesamples, const int &nLoops);
void call_kernel_PD_GPUv1_const(float *const d_input, float *const d_temp, unsigned char *const d_output_taps,
				const int &maxTaps, const int &nTimesamples, const float &signal_mean, const float &signal_sd);

} //namespace astroaccelerate
  
#endif
