#ifndef ASTRO_ACCELERATE_AA_DEVICE_MSD_OUTLIER_REJECTION_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_MSD_OUTLIER_REJECTION_KERNEL_HPP

namespace astroaccelerate {

void call_kernel_MSD_BLN_calculate_partials_2d_and_minmax_with_outlier_rejection(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, float *const d_MSD, const int &y_steps, const int &nTimesamples, const int &offset, const float &bln_sigma_constant);
void call_kernel_MSD_GPU_calculate_partials_2d_and_minmax(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, const int &y_steps, const int &nTimesamples, const int &offset);

void call_kernel_MSD_BLN_grid_outlier_rejection_new(const dim3 &grid_size, const dim3 &block_size, float *const d_input, float *const d_output, const int &size, const float &multiplier);
void call_kernel_MSD_BLN_grid_calculate_partials(const dim3 &grid_size, const dim3 &block_size, const int &threads, float const *const d_input, float *const d_output, const int &x_steps, const int &y_steps, const int &nColumns, const int &msd);
void call_kernel_MSD_BLN_grid_outlier_rejection(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, const int &size, const float &nElements, const float &multiplier);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_MSD_OUTLIER_REJECTION_KERNEL_HPP
