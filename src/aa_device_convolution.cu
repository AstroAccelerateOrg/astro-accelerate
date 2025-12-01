#include <vector>
#include "aa_device_convolution.hpp"

namespace astroaccelerate {

	void CONV_init(){
		//---------> Specific nVidia stuff
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
	}


	void forwardCustomFFT(float2 *d_filters, int FFT_size, int nFilters){
		dim3 gridSize(nFilters, 1, 1);
		dim3 blockSize(FFT_size/4, 1, 1);
		
		call_kernel_k_customFFT_GPU_forward(gridSize, blockSize, d_filters, d_filters, FFT_size);
	}


	void conv_OLS_customFFT(
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
	){
		dim3 gridSize(nConvolutions, 1, 1);
		dim3 blockSize(convolution_length/4, 1, 1);
		
		call_kernel_k_GPU_conv_OLS_via_customFFT(gridSize, blockSize, d_input_signal, d_output_plane, d_filters, signal_length, useful_part_size, offset, nConvolutions, nFilters, scale, convolution_length);
	}


	void convolve_signal_C2C(
			float2 *d_input_signal, 
			float *d_output_plane, 
			float2 *d_filters, 
			int64_t signal_length, 
			int64_t nFilters, 
			int64_t filter_halfwidth, 
			int64_t convolution_length, 
			float scale
	){
		int64_t useful_part_size = convolution_length - 2*filter_halfwidth + 1;
		int64_t nSegments        = (signal_length + useful_part_size - 1)/useful_part_size;
		conv_OLS_customFFT(d_input_signal, d_output_plane, d_filters, signal_length, convolution_length, useful_part_size, filter_halfwidth, nSegments, nFilters, scale);
	}

}
