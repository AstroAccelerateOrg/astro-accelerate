#ifndef ASTRO_ACCELERATE_AA_DEVICE_THRESHOLD_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_THRESHOLD_KERNEL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"

namespace astroaccelerate {

  /** Kernel wrapper function for THR_GPU_WARP kernel function. */
  void call_kernel_THR_GPU_WARP(const dim3 &grid_size, const dim3 &block_size,
				float const *const d_input, ushort *const d_input_taps, unsigned int *const d_output_list_DM,
				unsigned int *const d_output_list_TS, float *const d_output_list_SNR,
				unsigned int *const d_output_list_BW, int *const gmem_pos,
				const float &threshold, const int &nTimesamples, const int &offset, const int &shift, const int &max_list_size,
				const int &DIT_value);

  /** Kernel wrapper function for GPU_Threshold_for_periodicity_kernel_old kernel function. */
  void call_kernel_GPU_Threshold_for_periodicity_kernel_old(const dim3 &grid_size, const dim3 &block_size,
							    float const *const d_input, ushort *const d_input_harms,
							    float *const d_output_list, int *const gmem_pos, float *const d_MSD,
							    const float &threshold, const int &primary_size,
							    const int &secondary_size, const int &DM_shift,
							    const int &max_list_size, const int &DIT_value);
  /** Kernel wrapper function for GPU_Threshold_for_periodicity_kernel kernel function. */
  void call_kernel_GPU_Threshold_for_periodicity_kernel(const dim3 &grid_size, const dim3 &block_size,
							float const *const d_input, ushort *const d_input_harms,
							float *const d_output_list, int *const gmem_pos, float const *const d_MSD,
							const float &threshold, const int &primary_size, const int &secondary_size,
							const int &DM_shift, const int &max_list_size, const int &DIT_value);

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_DEVICE_THRESHOLD_KERNEL_HPP
