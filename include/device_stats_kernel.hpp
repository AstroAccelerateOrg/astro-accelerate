#ifndef ASTRO_ACCELERATE_DEVICE_STATS_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_STATS_KERNEL_HPP

void call_kernel_stats_kernel(dim3 block_size, dim3 grid_size, int smem_bytes, cudaStream_t stream,
			      int half_samps, float *d_sum, float *d_sum_square, float *d_signal_power);

#endif
