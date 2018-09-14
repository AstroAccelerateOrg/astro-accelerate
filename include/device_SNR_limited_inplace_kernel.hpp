#ifndef ASTRO_ACCELERATE_DEVICE_SNR_LIMITED_INPLACE_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_SNR_LIMITED_INPLACE_KERNEL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include "params.hpp"
#include "device_SPS_inplace_kernel.hpp"

void call_kernel_PD_ZC_GPU(float *d_input, float *d_output, int maxTaps, int nTimesamples, int nLoops);
void call_kernel_PD_GPUv1_const(float *d_input, float *d_temp, unsigned char *d_output_taps,
				int maxTaps, int nTimesamples, float signal_mean, float signal_sd);

#endif
