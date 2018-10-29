#ifndef ASTRO_ACCELERATE_DEVICE_ZERO_DM_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_ZERO_DM_KERNEL_HPP

void call_kernel_zero_dm_kernel(const dim3 &block_size, const dim3 &grid_size, const int &sharedMemory_size,
				unsigned short *const d_input, const int &nchans, const int &nsamp, const float &normalization_factor);

#endif
