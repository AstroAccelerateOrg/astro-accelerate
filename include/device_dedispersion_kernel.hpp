#ifndef ASTRO_ACCELERATE_DEDISPERSION_KERNEL_HPP
#define ASTRO_ACCELERATE_DEDISPERSION_KERNEL_HPP

#include "params.hpp"

//These device variables and definitions are needed by device_dedispersion_kernel.cu and device_load_data.cu
// Stores temporary shift values
#define ARRAYSIZE SDIVINT * SDIVINDM

__device__ __constant__ float dm_shifts[8192];
__device__ __constant__ int i_nsamp, i_nchans, i_t_processed_s;
__device__  __shared__ ushort2 f_line[UNROLLS][ARRAYSIZE + 2];

void call_kernel_shared_dedisperse_kernel(dim3 block_size, dim3 grid_size, int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep);
void call_kernel_shared_dedisperse_kernel_16(dim3 block_size, dim3 grid_size, int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep);
void call_kernel_cache_dedisperse_kernel(dim3 block_size, dim3 grid_size, int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep);

#endif

