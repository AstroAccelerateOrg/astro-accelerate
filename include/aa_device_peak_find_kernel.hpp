#ifndef ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_KERNEL_HPP

#include "aa_device_threshold_kernel.hpp"

namespace astroaccelerate {

void call_kernel_dilate_peak_find(const dim3 &grid_size, const dim3 &block_size,
				  float *const d_input, ushort *const d_input_taps,  unsigned int *const d_peak_list_DM,
				  unsigned int *const d_peak_list_TS, float *const d_peak_list_SNR, unsigned int *const d_peak_list_BW,
				  const int &width, const int &height, const int &offset, const float &threshold,
				  const int &max_peak_size, int *const gmem_pos, const int &shift, const int &DIT_value);

void call_kernel_dilate_peak_find_for_fdas(const dim3 &grid_size, const dim3 &block_size,
					   float *const d_input, float *const d_peak_list, float *const d_MSD, const int &width,
					   const int &height, const int &offset, const float &threshold,
					   const unsigned int &max_peak_size, unsigned int *const gmem_pos, const float &DM_trial);

void call_kernel_dilate_peak_find_for_periods_old(const dim3 &grid_size, const dim3 &block_size,
						  float *const d_input, ushort *const d_input_taps, float *const d_peak_list,
						  const int &width, const int &height, const int &offset,
						  const float &threshold,
						  const int &max_peak_size, int *const gmem_pos, float *const d_MDF,
						  const int &DM_shift, const int &DIT_value);

void call_kernel_dilate_peak_find_for_periods(const dim3 &grid_size, const dim3 &block_size,
					      float *const d_input, ushort *const d_input_taps,
					      float *const d_peak_list, const int &width,
					      const int &height, const int &offset, const float &threshold,
					      const int &max_peak_size,
					      int *const gmem_pos, float const *const d_MSD, const int &DM_shift, const int &DIT_value);

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_KERNEL_HPP
