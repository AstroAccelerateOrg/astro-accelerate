#ifndef ASTRO_ACCELERATE_DEVICE_HARMONIC_SUMMING_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_HARMONIC_SUMMING_KERNEL_HPP

void call_kernel_PHS_GPU_kernel_old(dim3 grid_size, dim3 block_size,
				    float const* d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics);

void call_kernel_PHS_GPU_kernel(dim3 grid_size, dim3 block_size,
				float const* d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics);

#endif
