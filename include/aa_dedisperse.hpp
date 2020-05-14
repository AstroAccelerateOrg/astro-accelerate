#ifndef ASTRO_ACCELERATE_AA_DEDISPERSE_HPP
#define ASTRO_ACCELERATE_AA_DEDISPERSE_HPP

#include <stdio.h>
#include <math.h>
#include <vector_types.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include "aa_params.hpp"
#include "aa_device_dedispersion_kernel.hpp"

namespace astroaccelerate {

  /**
   * \brief Function that performs the dedispersion on the GPU.
   * \brief Users should not need to interact with this function directly.
   */
  int dedisperse(int i, int t_processed, int *inBin, float *dmshifts, unsigned short *d_input, float *d_output, float *d_dm_shifts, int nchans, float *tsamp, float *dm_low, float *dm_step, int const*const ndms, int nbits, int failsafe);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEDISPERSE_HPP
