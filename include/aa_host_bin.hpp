#ifndef ASTRO_ACCELERATE_AA_HOST_BIN_HPP
#define ASTRO_ACCELERATE_AA_HOST_BIN_HPP

namespace astroaccelerate {

  /** \brief Performs binning on the host. */
  void bin(size_t binsize, float *bin_buffer,size_t inputsize, float *input_buffer, int nchans, int nsamps);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_HOST_BIN_HPP
