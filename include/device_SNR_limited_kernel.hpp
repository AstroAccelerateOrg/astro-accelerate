#ifndef ASTRO_ACCELERATE_DEVICE_SNR_LIMITED_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_SNR_LIMITED_KERNEL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include "params.hpp"

void call_kernel_SNR_GPU_limited(dim3 grid_size, dim3 block_size, float *d_FIR_input, float *d_SNR_output,
				 ushort *d_SNR_taps, float *d_MSD, int x_steps, int nTaps,
				 int nColumns, int offset);

#endif
