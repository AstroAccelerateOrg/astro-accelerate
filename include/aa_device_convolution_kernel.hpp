#ifndef ASTRO_ACCELERATE_AA_DEVICE_CONVOLUTION_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_CONVOLUTION_KERNEL_HPP

namespace astroaccelerate {

	void call_kernel_k_customFFT_GPU_forward(
			const dim3 &grid_size, 
			const dim3 &block_size, 
			float2 *d_input, 
			float2* d_output, 
			int FFT_size);

	void call_kernel_k_GPU_conv_OLS_via_customFFT(
			const dim3 &grid_size, 
			const dim3 &block_size, 
			float2 *d_input_signal, 
			float *d_output_plane, 
			float2 *d_filters, 
			int64_t signal_length, 
			int64_t useful_part_size, 
			int64_t offset, 
			int64_t nConvolutions, 
			int64_t nFilters,
			float scale,
			int64_t convolution_length);
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_CONVOLUTION_KERNEL_HPP

