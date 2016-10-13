#ifndef POWER_KERNEL_H_
#define POWER_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "AstroAccelerate/params.h"

<<<<<<< HEAD
//{{{ Set stretch
=======
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
__global__ void power_kernel(int half_samps, int acc, cufftComplex *d_signal_fft, float *d_signal_power)
{

	int t  = blockIdx.x * blockDim.x + threadIdx.x;
	
<<<<<<< HEAD
	if(t < half_samps) d_signal_power[t+acc*(half_samps)]= (d_signal_fft[t+1].x*d_signal_fft[t+1].x + d_signal_fft[t+1].y*d_signal_fft[t+1].y);
}

//}}}

=======
	if(t < half_samps) d_signal_power[t+acc*(half_samps)] = (d_signal_fft[t+1].x*d_signal_fft[t+1].x + d_signal_fft[t+1].y*d_signal_fft[t+1].y);
}

>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
#endif

