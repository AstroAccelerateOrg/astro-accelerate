#ifndef ASTRO_ACCELERATE_AA_DEVICE_STRETCH_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_STRETCH_KERNEL_HPP

namespace astroaccelerate {

  /** \brief Kernel wrapper function for stretch_kernel kernel function. */
  void call_kernel_stretch_kernel(const dim3 &block_size, const dim3 &grid_size, const int &smem_bytes, const cudaStream_t &stream,
				  const int &acc, const int &samps, const float &tsamp, float *const d_input, float *const d_output, const float &t_zero,
				  const float &multiplier, const float &tsamp_inverse);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_STRETCH_KERNEL_HPP
