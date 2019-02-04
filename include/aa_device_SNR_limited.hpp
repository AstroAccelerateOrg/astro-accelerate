#ifndef ASTRO_ACCELERATE_AA_DEVICE_SNR_LIMITED_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SNR_LIMITED_HPP

#include "aa_params.hpp"
#include "aa_device_SNR_limited_kernel.hpp"

namespace astroaccelerate {

  extern void SNR_limited_init(void);
  extern int SNR_limited(float *d_FIR_input, float *d_SNR_output, float *d_SNR_taps, float *d_MSD, int nTaps, int nDMs, int nTimesamples, int offset);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SNR_LIMITED_HPP

