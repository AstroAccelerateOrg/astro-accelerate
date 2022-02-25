#ifndef ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_HPP

namespace astroaccelerate {

extern int periodicity_simple_harmonic_summing(
  float *d_input,
  float *d_output_SNR,
  ushort *d_output_harmonics,
  float *d_MSD,
  int nTimesamples,
  int nDMs,
  int nHarmonics
);

extern int periodicity_greedy_harmonic_summing(
  float *d_input,
  float *d_output_SNR,
  ushort *d_output_harmonics,
  float *d_MSD,
  int nTimesamples,
  int nDMs,
  int nHarmonics,
  int enable_scalloping_loss_removal
);

extern int periodicity_presto_plus_harmonic_summing(
  float *d_input,
  float *d_output_SNR,
  ushort *d_output_harmonics,
  float *d_MSD,
  int nTimesamples,
  int nDMs,
  int nHarmonics,
  int enable_scalloping_loss_removal
);

extern int periodicity_presto_harmonic_summing(
  float *d_input,
  float *d_output_SNR,
  ushort *d_output_harmonics,
  float *d_MSD,
  int nTimesamples,
  int nDMs,
  int nHarmonics,
  int enable_scalloping_loss_removal
);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_HPP
