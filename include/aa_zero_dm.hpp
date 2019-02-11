#ifndef ASTRO_ACCELERATE_AA_ZERO_DM_HPP
#define ASTRO_ACCELERATE_AA_ZERO_DM_HPP

#include <time.h>
#include <math.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <vector_types.h>

#include "aa_params.hpp"
#include "aa_device_zero_dm_kernel.hpp"

namespace astroaccelerate {

  /**
   * \brief Function that performs zero_dm (without performing outlier rejection).
   * \details The user should not have to interact directly with this function.
   */
  void zero_dm(unsigned short *const d_input, const int nchans, const int nsamp, const int nbits);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_ZERO_DM_HPP
