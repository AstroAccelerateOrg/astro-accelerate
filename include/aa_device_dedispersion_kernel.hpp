#ifndef ASTRO_ACCELERATE_AA_DEDISPERSION_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEDISPERSION_KERNEL_HPP

#include "aa_params.hpp"
#include <vector_types.h>

namespace astroaccelerate {
  //These device variables and definitions are needed by device_dedispersion_kernel.cu and device_load_data.cu
  // Stores temporary shift values
#define ARRAYSIZE SDIVINT * SDIVINDM

  /** \brief Kernel wrapper function to set device constants for dedispersion_kernel kernel function. */
  void set_device_constants_dedispersion_kernel(const int &nchans, const int &length, const int &t_processed, const float *const dm_shifts);

  /** \brief Kernel wrapper function to set device constants for dedispersion_kernel kernel function. */
  void set_device_constants_dedispersion_kernel(const long int &length, const int &t_processed);

  /** \brief Kernel wrapper function for shared_dedisperse_kernel kernel function. */
  void call_kernel_shared_dedisperse_kernel(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);

  /** \brief Kernel wrapper function for shared_dedisperse_kernel_16 kernel function. */
  void call_kernel_shared_dedisperse_kernel_16(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);

  /** \brief Kernel wrapper function for cache_dedisperse_kernel kernel function. */
  void call_kernel_cache_dedisperse_kernel(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_DEDISPERSION_KERNEL_HPP

