#ifndef ASTRO_ACCELERATE_DEVICE_ZERO_DM_OUTLIERS_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_ZERO_DM_OUTLIERS_KERNEL_HPP

void call_kernel_zero_dm_outliers_kernel_one(const dim3 &          block_size,
                                             const dim3 &          grid_size,
                                             unsigned short *const d_input,
                                             const int &           nchans,
                                             const int &           nsamp);

void call_kernel_zero_dm_outliers_kernel_two(const dim3 &          block_size,
                                             const dim3 &          grid_size,
                                             unsigned short *const d_input,
                                             const int &           nchans,
                                             const int &           nsamp);

#endif
