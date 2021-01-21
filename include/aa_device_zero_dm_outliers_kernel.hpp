#ifndef ASTRO_ACCELERATE_AA_DEVICE_ZERO_DM_OUTLIERS_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_ZERO_DM_OUTLIERS_KERNEL_HPP

namespace astroaccelerate {

void call_kernel_zero_dm_outliers_kernel_channels(
	const dim3 &grid_size, 
	const dim3 &block_size, 
	const int &smem, 
	const cudaStream_t &stream, 
	unsigned short *const d_input, 
	const int &nchans, 
	const float &outlier_sigma, 
	const int &nbits,
	float *normalization_factor);
	
  /** \brief Kernel wrapper function for zero_dm_outliers_kernel_one kernel function. */
  void call_kernel_zero_dm_outliers_kernel_one(const dim3 &block_size, const dim3 &grid_size,
					       unsigned short *const d_input, const int &nchans, const int &nsamp);

  /** \brief Kernel wrapper function for zero_dm_outliers_kernel_two kernel function. */
  void call_kernel_zero_dm_outliers_kernel_two(const dim3 &block_size, const dim3 &grid_size,
					       unsigned short *const d_input, const int &nchans, const int &nsamp);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_ZERO_DM_OUTLIERS_KERNEL_HPP
