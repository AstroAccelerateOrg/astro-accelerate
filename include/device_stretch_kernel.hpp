#ifndef ASTRO_ACCELERATE_DEVICE_STRETCH_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_STRETCH_KERNEL_HPP

void call_kernel_stretch_kernel(dim3 block_size, dim3 grid_size, int smem_bytes, cudaStream_t stream,
				int acc, int samps, float tsamp, float *d_input, float *d_output, float t_zero,
				float multiplier, float tsamp_inverse);

#endif
