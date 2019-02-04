#ifndef ASTRO_ACCELERATE_AA_DEVICE_PERIODS_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_PERIODS_HPP

namespace astroaccelerate {

extern void GPU_periodicity(int range, int nsamp, int max_ndms, int processed, float sigma_cutoff, float ***output_buffer, int const*const ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp, int nHarmonics, bool candidate_algorithm, bool enable_msd_baseline_noise, float bln_sigma_constant);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_PERIODS_HPP








