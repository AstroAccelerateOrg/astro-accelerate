#ifndef ASTRO_ACCELERATE_AA_HOST_ACCELERATION_HPP
#define ASTRO_ACCELERATE_AA_HOST_ACCELERATION_HPP

namespace astroaccelerate {

  /**
   * \brief Function that performs accelerated search.
   * \brief The user should not interact with this function directly and instead use fdas.
   */
  void acceleration(int range, int nsamp, int max_ndms, int processed, int num_boots, int num_trial_bins, int navdms, float narrow, float wide, int nsearch, float aggression, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_HOST_ACCELERATION_HPP

