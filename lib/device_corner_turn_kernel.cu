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

	//d_output[(size_t)(c * nsamp) + (size_t)t] = (float) __ldg(&d_input[(size_t)(t * nchans) + (size_t)c]);
	d_output[(size_t)(c * nsamp) + (size_t)t] = (float) (d_input[(size_t)(t * nchans) + (size_t)c]);

}

__global__ void swap(unsigned short *d_input, float *d_output, int nchans, int nsamp)
{

	int t = blockIdx.x * blockDim.x + threadIdx.x;
	int c = blockIdx.y * blockDim.y + threadIdx.y;

	//d_input[(size_t)(c * nsamp) + (size_t)t] = (unsigned short) __ldg(&d_output[(size_t)(c * nsamp) + (size_t)t]);
	d_input[(size_t)(c * nsamp) + (size_t)t] = (unsigned short) (d_output[(size_t)(c * nsamp) + (size_t)t]);

}

//}}}

#endif

