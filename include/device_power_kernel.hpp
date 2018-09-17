#ifndef ASTRO_ACCELERATE_DEVICE_POWER_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_POWER_KERNEL_HPP

void call_kernel_power_kernel(dim3 block_size, dim3 grid_size, int smem_bytes, cudaStream_t stream,
			      int half_samps, int acc, cufftComplex *d_signal_fft, float *d_signal_power);
void call_kernel_GPU_simple_power_and_interbin_kernel(dim3 grid_size, dim3 block_size,
						      float2 *d_input_complex, float *d_output_power, float *d_output_interbinning, int nTimesamples, float norm);

#endif
