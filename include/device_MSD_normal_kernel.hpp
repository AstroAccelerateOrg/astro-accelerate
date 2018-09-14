#ifndef ASTRO_ACCELERATE_DEVICE_MSD_NORMAL_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_MSD_NORMAL_KERNEL_HPP

void call_kernel_MSD_GPU_limited(dim3 grid_size, dim3 block_size, float const* d_input, float* d_output, int y_steps, int nTimesamples, int offset);

#endif
