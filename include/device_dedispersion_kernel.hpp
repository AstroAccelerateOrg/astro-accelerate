#ifndef ASTRO_ACCELERATE_DEDISPERSION_KERNEL_HPP
#define ASTRO_ACCELERATE_DEDISPERSION_KERNEL_HPP

#include "params.hpp"
#include <vector_types.h>
//These device variables and definitions are needed by device_dedispersion_kernel.cu and device_load_data.cu
// Stores temporary shift values
#define ARRAYSIZE SDIVINT * SDIVINDM

void set_device_constants_dedispersion_kernel(const int &nchans, const int &length, const int &t_processed, const float *const dm_shifts);
void set_device_constants_dedispersion_kernel(const long int &length, const int &t_processed);
void call_kernel_shared_dedisperse_kernel(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);
void call_kernel_shared_dedisperse_kernel_16(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);
void call_kernel_cache_dedisperse_kernel(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);

#endif

