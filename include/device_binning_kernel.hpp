#ifndef ASTRO_ACCELERATE_DEVICE_BINNING_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_BINNING_KERNEL_HPP

void call_kernel_bin(dim3 num_blocks, dim3 threads_per_block, unsigned short *d_input, float *d_output, int in_nsamp);
void call_kernel_DiT_GPU_v2(dim3 gridSize, dim3 blockSize, float const* d_input, float *d_output, unsigned int nDMs, unsigned int nTimesamples, unsigned int dts);

#endif
