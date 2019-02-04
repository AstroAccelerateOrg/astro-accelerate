#ifndef ASTRO_ACCELERATE_AA_DEVICE_RFI_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_RFI_KERNEL_HPP

namespace astroaccelerate {

  /** \brief Kernel wrapper function for call to rfi_gpu_kernel kernel function. */
  void call_kernel_rfi_gpu_kernel(const dim3 &grid_size, const dim3 &block_size,
				  unsigned short *const d_input,
				  const int &nchans, const int &nsamp);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_RFI_KERNEL_HPP
