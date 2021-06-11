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
    float *normalization_factor
  );

  void call_kernel_post_DDTR_normalization_dm(
      dim3 &nBlocks_per_grid, 
      dim3 &nThreads_per_block, 
      int shared_memory, 
      float *d_input, 
      size_t nTimesamples, 
      size_t nDMs,
      float normalization_factor,
      int nbits
  );
  
  void call_kernel_zero_dm_normalization_dm(
      dim3 &nBlocks_per_grid, 
      dim3 &nThreads_per_block, 
      int shared_memory_size, 
      unsigned short *d_input, 
      size_t nTimesamples, 
      size_t nDMs,
      float normalization_factor,
      int nbits
  );

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_ZERO_DM_KERNEL_HPP
