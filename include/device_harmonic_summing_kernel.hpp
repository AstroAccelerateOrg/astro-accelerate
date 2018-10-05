#ifndef ASTRO_ACCELERATE_DEVICE_HARMONIC_SUMMING_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_HARMONIC_SUMMING_KERNEL_HPP

void call_kernel_PHS_GPU_kernel_old(const dim3 &grid_size, const dim3 &block_size,
				    float const *const d_input, float *const d_output_SNR, ushort *const d_output_harmonics, float *const d_MSD, const int &nTimesamples, const int &nSpectra, const int &nHarmonics);

void call_kernel_PHS_GPU_kernel(const dim3 &grid_size, const dim3 &block_size,
				float const *const d_input, float *const d_output_SNR, ushort *const d_output_harmonics, float *const d_MSD, const int &nTimesamples, const int &nSpectra, const int &nHarmonics);

#endif
