#ifndef ASTRO_ACCELERATE_AA_DEVICE_SNR_LIMITED_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SNR_LIMITED_KERNEL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"

namespace astroaccelerate {

  /** \brief Kernel wrapper function for SNR_GPU_limited kernel function. */
  void call_kernel_SNR_GPU_limited(const dim3 &grid_size, const dim3 &block_size, float *const d_FIR_input, float *const d_SNR_output,
				   ushort *const d_SNR_taps, float *const d_MSD, const int &x_steps, const int &nTaps,
				   const int &nColumns, const int &offset);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SNR_LIMITED_KERNEL_HPP
