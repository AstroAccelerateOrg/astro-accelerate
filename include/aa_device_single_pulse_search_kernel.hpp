#ifndef ASTRO_ACCELERATE_AA_DEVICE_SINGLE_PULSE_SEARCH_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SINGLE_PULSE_SEARCH_KERNEL_HPP

#include "aa_params.hpp"

namespace astroaccelerate {

  /** \brief Kernel wrapper function for PD_SEARCH_GPU kernel function. */
  void call_kernel_PD_SEARCH_GPU(const dim3 &grid_size, const dim3 &block_size, const int &sm_size,
				 float const *const d_input, float *const d_output, float *const d_output_taps, float *const d_MSD, const int &maxTaps, const int &nTimesamples);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SINGLE_PULSE_SEARCH_KERNEL_HPP

