#ifndef ASTRO_ACCELERATE_AA_DEVICE_CONVOLUTION_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_CONVOLUTION_HPP

#include "aa_device_convolution_kernel.hpp"

namespace astroaccelerate {
	extern void CONV_init(void);

	extern void forwardCustomFFT(float2 *d_filters, int FFT_size, int nFilters);

	extern void conv_OLS_customFFT(
			float2 *d_input_signal, 
			float *d_output_plane, 
			float2 *d_filters, 
			int64_t signal_length, 
			int64_t convolution_length, 
			int64_t useful_part_size, 
			int64_t offset, 
			int64_t nConvolutions, 
			int64_t nFilters, 
			float scale
	);

	extern void convolve_signal_C2C(
			float2 *d_input_signal, 
			float *d_output_plane, 
			float2 *d_filters, 
			int64_t signal_length, 
			int64_t nFilters, 
			int64_t filter_halfwidth, 
			int64_t convolution_length, 
			float scale
	);
}
#endif
