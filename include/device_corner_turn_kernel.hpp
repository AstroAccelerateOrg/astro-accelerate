#ifndef ASTRO_ACCELERATE_DEVICE_CORNER_TURN_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_CORNER_TURN_KERNEL_HPP

void call_kernel_simple_corner_turn_kernel(dim3 block_size, dim3 grid_size, unsigned short *d_input, float *d_output, int primary_size, int secondary_size);
void call_kernel_simple_corner_turn_kernel(dim3 block_size, dim3 grid_size, float *d_input, float *d_output, int primary_size, int secondary_size);
void call_kernel_corner_turn_SM_kernel(dim3 grid_size, dim3 block_size, float const* d_input, float *d_output, int primary_size, int secondary_size);
void call_kernel_swap(dim3 block_size, dim3 grid_size, unsigned short *d_input, float *d_output, int nchans, int nsamp);
#endif
