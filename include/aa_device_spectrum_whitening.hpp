#ifndef ASTRO_ACCELERATE_AA_DEVICE_SPECTRUM_WHITENING_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SPECTRUM_WHITENING_HPP

#include <cufft.h>

namespace astroaccelerate {

extern void spectrum_whitening_SGP1(float *d_input, unsigned long int nSamples, int nDMs, cudaStream_t &stream);


} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SPECTRUM_WHITENING_HPP

