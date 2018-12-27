#ifndef ASTRO_ACCELERATE_DEDISPERSE_HPP
#define ASTRO_ACCELERATE_DEDISPERSE_HPP

#include <stdio.h>
#include <math.h>
#include <vector_types.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include "params.hpp"
#include "device_dedispersion_kernel.hpp"

namespace astroaccelerate {

void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, unsigned short *d_input, float *d_output, int nchans, float *tsamp, float *dm_low, float *dm_step, int const*const ndms, int nbits, int failsafe);

} //namespace astroaccelerate
  
#endif /* ASTRO_ACCELERATE_DEDISPERSE_HPP */
