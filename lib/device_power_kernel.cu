#ifndef POWER_KERNEL_H_
#define POWER_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "headers/params.h"

//{{{ Set stretch
__global__ void power_kernel(int half_samps, int acc, cufftComplex *d_signal_fft, float *d_signal_power)
{

	int t = blockIdx.x * blockDim.x + threadIdx.x;

	if (t < half_samps)
		d_signal_power[t + acc * ( half_samps )] = ( d_signal_fft[t + 1].x * d_signal_fft[t + 1].x + d_signal_fft[t + 1].y * d_signal_fft[t + 1].y );
}

//}}}


__global__ void GPU_simple_power_and_interbin_kernel(float2 *d_input_complex, float *d_output_power, float *d_output_interbinning, int nTimesamples){
	int half  = (nTimesamples>>1);
	int pos_x = blockIdx.x*blockDim.x;
	int pos_y = blockIdx.y*nTimesamples;
	
	float2 A, B;
	if ( (pos_x+threadIdx.x) < half) {
		A = d_input_complex[pos_y + pos_x + threadIdx.x];
		B = d_input_complex[pos_y + pos_x + threadIdx.x + 1];
		
		d_output_power[(pos_y>>1) + pos_x + threadIdx.x] = A.x*A.x + A.y*A.y;
		d_output_interbinning[pos_y + 2*(pos_x + threadIdx.x)] = A.x*A.x + A.y*A.y;
		d_output_interbinning[pos_y + 2*(pos_x + threadIdx.x) + 1] = 0.616850275f*( (A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y) );
	}
}

#endif

