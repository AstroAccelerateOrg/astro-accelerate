#ifndef ASTRO_ACCELERATE_DEVICE_SET_STRETCH_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_SET_STRETCH_KERNEL_HPP

void call_kernel_set_stretch_kernel(const dim3 &block_size, const dim3 &grid_size,
				    const int &smem_bytes, const cudaStream_t &stream,
				    const int &samps, const float &mean, float *const d_input);

#endif
