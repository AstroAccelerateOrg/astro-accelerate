#ifndef ASTRO_ACCELERATE_AA_DEVICE_CONVOLUTION_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_CONVOLUTION_HPP

#include "aa_device_convolution_kernel.hpp"

namespace astroaccelerate {
	extern void CONV_init(void);

	extern void forwardCustomFFT(float2 *d_filters, int FFT_size, int nFilters);

	extern void conv_OLS_customFFT(float2 *d_input_signal, float *d_output_plane, float2 *d_filters, int signal_length, int convolution_length, int useful_part_size, int offset, int nConvolutions, int nFilters, float scale);

	extern void convolve_signal_C2C(float2 *d_input_signal, float *d_output_plane, float2 *d_filters, int signal_length, int nFilters, int filter_halfwidth, int convolution_length, float scale);
}
#endif
