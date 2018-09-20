#ifndef ASTRO_ACCELERATE_DEVICE_MSD_OUTLIER_REJECTION_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_MSD_OUTLIER_REJECTION_KERNEL_HPP

void call_kernel_MSD_BLN_pw_rejection_normal(dim3 grid_size, dim3 block_size, float const* d_input, float *d_output, float *d_MSD, int y_steps, int nTimesamples, int offset, float bln_sigma_constant);
void call_kernel_MSD_BLN_grid_outlier_rejection_new(dim3 grid_size, dim3 block_size, float *d_input, float *d_output, int size, float multiplier);
void call_kernel_MSD_BLN_grid_calculate_partials(dim3 grid_size, dim3 block_size, int threads, float const* d_input, float *d_output, int x_steps, int y_steps, int nColumns, int msd);
void call_kernel_MSD_BLN_grid_outlier_rejection(dim3 grid_size, dim3 block_size, float const* d_input, float *d_output, int size, float nElements, float multiplier);

#endif
