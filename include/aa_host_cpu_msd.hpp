#ifndef ASTRO_ACCELERATE_AA_HOST_CPU_MSD_HPP
#define ASTRO_ACCELERATE_AA_HOST_CPU_MSD_HPP

#include <stdlib.h>

namespace astroaccelerate {

  /** \brief Performs mean and standard deviation on CPU. */
  void call_cpu_msd(float* h_new_bandpass, unsigned short const*const input_buffer, size_t nchans, size_t nsamples);

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_HOST_CPU_MSD_HPP

