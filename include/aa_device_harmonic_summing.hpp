#ifndef ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_HPP

namespace astroaccelerate {

extern void periodicity_simple_harmonic_summing_old(float *d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics);

extern void periodicity_simple_harmonic_summing(float *d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_HPP
