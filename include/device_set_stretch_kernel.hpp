#ifndef ASTRO_ACCELERATE_DEVICE_SET_STRETCH_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_SET_STRETCH_KERNEL_HPP

void call_kernel_set_stretch_kernel(dim3 block_size, dim3 grid_size,
				    int smem_bytes, cudaStream_t stream,
				    int samps, float mean, float *d_input);

#endif
