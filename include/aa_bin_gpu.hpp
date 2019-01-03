#ifndef ASTRO_ACCELERATE_AA_BIN_GPU_HPP
#define ASTRO_ACCELERATE_AA_BIN_GPU_HPP

#include <stdio.h>
#include <vector_types.h>

#include "aa_params.hpp"
#include "aa_device_binning_kernel.hpp"
#include "aa_device_corner_turn_kernel.hpp"

namespace astroaccelerate {

  /**
   * Functions to perform binning.
   */

  void bin_gpu(unsigned short *const d_input, float *const d_output, const int nchans, const int nsamp);
  int GPU_DiT_v2_wrapper(float *d_input, float *d_output, int nDMs, int nTimesamples);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_BIN_GPU_HPP
