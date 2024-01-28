#ifndef ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_KERNEL_HPP

#include "aa_device_threshold_kernel.hpp"

namespace astroaccelerate {

  /** \brief Kernel wrapper function for dilate_peak_find kernel function. */
  void call_kernel_dilate_peak_find(const dim3 &grid_size, const dim3 &block_size,
				    float *const d_input, ushort *const d_input_taps,  unsigned int *const d_peak_list_DM,
				    unsigned int *const d_peak_list_TS, float *const d_peak_list_SNR, unsigned int *const d_peak_list_BW,
				    const int &width, const int &height, const int &offset, const float &threshold,
				    const int &max_peak_size, int *const gmem_pos, const int &shift, const int &DIT_value);

  /** \brief Kernel wrapper function for dilate_peak_find_for_fdas kernel function. */
  void call_kernel_dilate_peak_find_for_fdas(const dim3 &grid_size, const dim3 &block_size,
					     float *const d_input, float *const d_peak_list, float *const d_MSD, const int &width,
					     const int &height, const int &offset, const float &threshold,
					     const unsigned int &max_peak_size, unsigned int *const gmem_pos, const float &DM_trial);
						 
	void call_kernel_dilate_peak_find_for_fdas_harm(
		const dim3 &gridSize, 
		const dim3 &blockDim, 
		float *d_peak_list, 
		float *d_ffdot_max, 
		float *d_ffdot_SNR, 
		ushort *d_ffdot_harm, 
		size_t nFreq, 
		size_t nAcc, 
		int half_plane, 
		float threshold, 
		unsigned int max_peak_size, 
		unsigned int *const gmem_peak_pos, 
		float DM_trial
	);

  /** \brief Kernel wrapper function for peak_find_for_periods_old kernel function. */
  void call_kernel_peak_find_for_periodicity_normal(
      const dim3 &grid_size,
      const dim3 &block_size,
      float *const d_input,
      ushort *const d_input_taps,
      float *const d_peak_list,
      const int &nTimesamples,
      const int &nDMs,
      const int &offset,
      const float &threshold,
      const int &max_peak_size,
      int *const gmem_pos,
      float const *const d_MSD,
      const int &DM_shift,
      const int &DIT_value
  );

  /** \brief Kernel wrapper function for peak_find_for_periods kernel function. */
  void call_kernel_peak_find_for_periodicity_transposed(
      const dim3 &grid_size,
      const dim3 &block_size,
      float *const d_input,
      ushort *const d_input_taps,
      float *const d_peak_list,
      const int &width,
      const int &height,
      const int &offset,
      const float &threshold,
      const int &max_peak_size,
      int *const gmem_pos,
      float const *const d_MSD,
      const int &DM_shift,
      const int &DIT_value
  );

  void call_kernel_peak_find_list(const dim3 &grid_size, const dim3 &block_size, float *const d_input, const int width, const int height, const float &threshold, int *const gmem_pos, const int &shift, const int &DIT_value, ushort *const d_input_taps, const int &max_peak_size, unsigned int *const d_peak_list_DM, unsigned int *const d_peak_list_TS, float *const d_peak_list_SNR, unsigned int *const d_peak_list_BW);

  void call_gpu_Filter_peaks(unsigned int *new_peak_list_DM, unsigned int *new_peak_list_TS, unsigned int *new_peak_list_BW, float *new_peak_list_SNR, unsigned int *d_peak_list_DM, unsigned int *d_peak_list_TS, unsigned int *d_peak_list_BW, float *d_peak_list_SNR, unsigned int nElements, unsigned int max_distance, int max_list_pos, int *gmem_pos);

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_KERNEL_HPP
