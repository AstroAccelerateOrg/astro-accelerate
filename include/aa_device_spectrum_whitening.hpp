#ifndef ASTRO_ACCELERATE_AA_DEVICE_SPECTRUM_WHITENING_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SPECTRUM_WHITENING_HPP

#include <stdio.h>
#include <vector>
#include <cufft.h>

namespace astroaccelerate {

extern void create_dered_segment_sizes_prefix_sum(std::vector<int> *segment_sizes, int min_segment_length, int max_segment_length, size_t nSamples);

extern void spectrum_whitening_SGP1(float *d_input, unsigned long int nSamples, int nDMs, cudaStream_t &stream);

extern void spectrum_whitening_SGP2(float2 *d_input, size_t nSamples, int nDMs, bool enable_median, cudaStream_t &stream);


} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SPECTRUM_WHITENING_HPP

