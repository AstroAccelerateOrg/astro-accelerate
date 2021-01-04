#ifndef ASTRO_ACCELERATE_AA_DEVICE_ZERO_DM_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_ZERO_DM_KERNEL_HPP

namespace astroaccelerate {

  /** \brief Kernel wrapper function for zero_dm_kernel kernel function. */
  void call_kernel_zero_dm_kernel(
    const dim3 &block_size, 
    const dim3 &grid_size, 
    const int &sharedMemory_size,
    unsigned short *const d_input, 
    const int &nchans, 
    const int &nsamp, 
    const int &nbits,
    const float &normalization_factor);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_ZERO_DM_KERNEL_HPP
