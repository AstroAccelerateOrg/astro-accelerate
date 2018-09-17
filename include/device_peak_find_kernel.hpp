#ifndef ASTRO_ACCELERATE_DEVICE_PEAK_FIND_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_PEAK_FIND_KERNEL_HPP

#include "device_threshold_kernel.hpp"


void call_kernel_dilate_peak_find(dim3 grid_size, dim3 block_size,
				  const float *d_input, ushort* d_input_taps,  unsigned int *d_peak_list_DM,
				  unsigned int *d_peak_list_TS, float *d_peak_list_SNR, unsigned int *d_peak_list_BW,
				  const int width, const int height, const int offset, const float threshold,
				  int max_peak_size, int *gmem_pos, int shift, int DIT_value);

void call_kernel_dilate_peak_find_for_fdas(dim3 grid_size, dim3 block_size,
					   const float *d_input, float *d_peak_list, float *d_MSD, const int width,
					   const int height, const int offset, const float threshold,
					   unsigned int max_peak_size, unsigned int *gmem_pos, float DM_trial);

void call_kernel_dilate_peak_find_for_periods_old(dim3 grid_size, dim3 block_size,
						  const float *d_input, ushort* d_input_taps, float *d_peak_list,
						  const int width, const int height, const int offset,
						  const float threshold,
						  int max_peak_size, int *gmem_pos, float *d_MDF,
						  int DM_shift, int DIT_value);

void call_kernel_dilate_peak_find_for_periods(dim3 grid_size, dim3 block_size,
					      const float *d_input, ushort* d_input_taps,
					      float *d_peak_list, const int width,
					      const int height, const int offset, const float threshold,
					      int max_peak_size,
					      int *gmem_pos, float const* d_MSD, int DM_shift, int DIT_value);
#endif
