#ifndef ASTRO_ACCELERATE_AA_HOST_STATISTICS_HPP
#define ASTRO_ACCELERATE_AA_HOST_STATISTICS_HPP

namespace astroaccelerate {

  /** \brief Function that prints statistics as provided to it. */
  void statistics(char *string, int i, cudaStream_t stream, double *in_time, double *out_time, int maxshift, int total_ndms, int nchans, int nsamp, float tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_HOST_STATISTICS_HPP
