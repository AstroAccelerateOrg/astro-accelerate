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

#endif

