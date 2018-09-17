#ifndef ASTRO_ACCELERATE_DEVICE_SINGLE_PULSE_SEARCH_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_SINGLE_PULSE_SEARCH_KERNEL_HPP

#include "params.hpp"

void call_kernel_PD_SEARCH_GPU(dim3 grid_size, dim3 block_size, int sm_size,
			       float const* d_input, float *d_output, float *d_output_taps, float *d_MSD, int maxTaps, int nTimesamples);

#endif

