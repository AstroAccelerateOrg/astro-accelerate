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

  /** \brief Kernel wrapper function to set device constants for dedispersion_kernel kernel function. */
  void set_device_constants_dedispersion_kernel(const int &nchans, const long int &length, const int &t_processed);

  /** \brief Kernel wrapper function for shared_dedisperse_kernel kernel function. */
  void call_kernel_shared_dedisperse_kernel(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);
  
	/** \brief Kernel wrapper function for dedispersion GPU kernel which works with number of channels greater then 8192. */
	void call_kernel_shared_dedisperse_kernel_nchan8192p(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, float *const d_dm_shifts, const float &mstartdm, const float &mdmstep);

  /** \brief Kernel wrapper function for shared_dedisperse_kernel_16 kernel function. */
  void call_kernel_shared_dedisperse_kernel_16(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);
  
	/** \brief Kernel wrapper function for dedispersion kernel which works with 16bit data and when number of channels is greater than 8192 kernel function. */
	void call_kernel_shared_dedisperse_kernel_16_nchan8192p(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, float *const d_dm_shifts, const float &mstartdm, const float &mdmstep);

  /** \brief Kernel wrapper function for cache_dedisperse_kernel kernel function. */
  void call_kernel_cache_dedisperse_kernel(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);
  
	void call_kernel_cache_dedisperse_kernel_nchan8192p(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, float *const d_dm_shifts, const float &mstartdm, const float &mdmstep);

	/** \brief kernel wrapper for dedispersion kernel which works with 4bit data in case when the number of channels is less than 4096*/
	void call_kernel_shared_dedisperse_kernel_4bit(const dim3 &block_size, const dim3 &grid_size, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep);

        /** \brief Kernel wrapper for dedispersion kernel which works with 4bit data in case when the nchans is greater than 4096*/
        void call_kernel_shared_dedisperse_kernel_4bit_4096chan(const dim3 &block_size, const dim3 &grid_size, unsigned short *const d_input, float *const d_output, float *const d_dm_shifts, const float &mstartdm, const float &mdmstep, const int nchans);

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_DEDISPERSION_KERNEL_HPP

