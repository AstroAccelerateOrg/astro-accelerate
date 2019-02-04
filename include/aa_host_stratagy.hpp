#ifndef ASTRO_ACCELERATE_AA_HOST_STRATAGY_HPP
#define ASTRO_ACCELERATE_AA_HOST_STRATAGY_HPP

namespace astroaccelerate {

  /**
   * \warning This function does not look like it has an implementation since there is no corresponding file called aa_host_stratagy.cpp or similar.
   */
  void stratagy(int *maxshift, int *max_samps, int *num_tchunks, int *max_ndms, int *total_ndms, float *max_dm, float power, int nchans, int nsamp, float fch1, float foff, float tsamp, int range, float *user_dm_low, float *user_dm_high, float *user_dm_step, float **dm_low, float **dm_high, float **dm_step, int **ndms, float **dmshifts, int *inBin, int ***t_processed, size_t *gpu_memory, int enable_analysis);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_HOST_STRATAGY_HPP
