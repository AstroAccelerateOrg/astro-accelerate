#ifndef ASTRO_ACCELERATE_DEVICE_ZERO_DM_OUTLIERS_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_ZERO_DM_OUTLIERS_KERNEL_HPP

void call_kernel_zero_dm_outliers_kernel_one(dim3 block_size, dim3 grid_size,
					     unsigned short *d_input, int nchans, int nsamp);

void call_kernel_zero_dm_outliers_kernel_two(dim3 block_size, dim3 grid_size,
					     unsigned short *d_input, int nchans, int nsamp);

#endif
