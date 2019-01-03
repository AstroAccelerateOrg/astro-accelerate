#ifndef ASTRO_ACCELERATE_AA_DEVICE_BINNING_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_BINNING_KERNEL_HPP

namespace astroaccelerate {

  /** \brief Kernel wrapper function for bin kernel function. */
  void call_kernel_bin(const dim3 &num_blocks, const dim3 &threads_per_block, unsigned short *const d_input, float *const d_output, const int &in_nsamp);

  /** \brief Kernel wrapper function for DiT_GPU_v2 kernel function. */
  void call_kernel_DiT_GPU_v2(const dim3 &gridSize, const dim3 &blockSize, float const *const d_input, float *const d_output, const unsigned int &nDMs, const unsigned int &nTimesamples, const unsigned int &dts);
  
} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_BINNING_KERNEL_HPP
