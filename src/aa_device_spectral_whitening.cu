#include <stdio.h>
#include <cufft.h>
#include "aa_params.hpp"
#include "aa_device_spectrum_whitening_kernel.hpp"

// define for debug info
//#define POWER_DEBUG

//{{{ Dopler Stretch 

namespace astroaccelerate {

	void spectrum_whitening_SGP1(float *d_input, unsigned long int nSamples, int nDMs, cudaStream_t &stream){
		dim3 grid_size((nSamples + 128 - 1)/128, nDMs , 1);
		dim3 block_size(128, 1, 1);
		call_kernel_spectrum_whitening_SGP1(
			block_size, 
			grid_size, 
			0, 
			stream, 
			d_input, 
			nSamples
		);
	}

} //namespace astroaccelerate
