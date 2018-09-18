#ifndef ASTRO_ACCELERATE_DEVICE_ZERO_DM_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_ZERO_DM_KERNEL_HPP

void call_kernel_zero_dm_kernel(dim3 block_size, dim3 grid_size,
				unsigned short *d_input, int nchans, int nsamp, float normalization_factor);

#endif
