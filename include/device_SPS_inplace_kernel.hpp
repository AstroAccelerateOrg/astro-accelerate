#ifndef ASTROACCELERATE_DEVICE_SPS_INPLACE_KERNEL_HPP
#define ASTROACCELERATE_DEVICE_SPS_INPLACE_KERNEL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include "params.hpp"

//This device constant is needed by the SPS module and the SNR module
__device__ __constant__ float c_sqrt_taps[PD_MAXTAPS + 1];

void call_kernel_PD_ZC_GPU_KERNEL(dim3 grid_size, dim3 block_size, float *d_input, float *d_output,
				  int maxTaps, int nTimesamples, int nLoops);
void call_kernel_PD_INPLACE_GPU_KERNEL(dim3 grid_size, dim3 block_size, int SM_size, float *d_input,
				       float *d_temp, unsigned char *d_output_taps, float *d_MSD,
				       int maxTaps, int nTimesamples);

#endif
