#ifndef ASTRO_ACCELERATE_AA_DEVICE_SPS_INPLACE_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SPS_INPLACE_HPP

#include "aa_params.hpp"
#include "aa_device_SPS_inplace_kernel.hpp"

namespace astroaccelerate {

  extern void PD_SEARCH_INPLACE_init(void);
  extern int PD_SEARCH_INPLACE(float *d_input, unsigned char *d_output_taps, float *d_MSD, int maxTaps, int nDMs, int nTimesamples);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SPS_INPLACE_HPP

