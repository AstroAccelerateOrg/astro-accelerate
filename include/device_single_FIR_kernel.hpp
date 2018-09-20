#ifndef ASTRO_ACCELERATE_DEVICE_SINGLE_FIR_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_SINGLE_FIR_KERNEL_HPP

void call_kernel_PD_FIR_GPU(dim3 grid_size, dim3 block_size, int SM_size, float const* d_input, float *d_output, int nTaps, int nLoops, int nTimesamples);
void call_kernel_PD_FIR_GPUv1(dim3 grid_size, dim3 block_size, int SM_size, float const* d_input, float *d_output, int nTaps, int nLoops, unsigned int nTimesamples);
void call_kernel_Fir_L1(dim3 grid_size, dim3 block_size, float const* d_input, float* d_output, int nTaps, int nTimesamples);

#endif
