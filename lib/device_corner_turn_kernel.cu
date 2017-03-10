#ifndef CORNERTURN_KERNEL_H_
#define CORNERTURN_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

//{{{ corner_turn
__global__ void simple_corner_turn_kernel(unsigned short *d_input, float *d_output, int nchans, int nsamp)
{

	int t = blockIdx.x * blockDim.x + threadIdx.x;
	int c = blockIdx.y * blockDim.y + threadIdx.y;

	d_output[c * nsamp + t] = (float) __ldg(&d_input[t * nchans + c]);

}

__global__ void swap(unsigned short *d_input, float *d_output, int nchans, int nsamp)
{

	int t = blockIdx.x * blockDim.x + threadIdx.x;
	int c = blockIdx.y * blockDim.y + threadIdx.y;

	d_input[c * nsamp + t] = (unsigned short) __ldg(&d_output[c * nsamp + t]);

}

//}}}

#endif

