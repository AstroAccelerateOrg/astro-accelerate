#ifndef SET_STRETCH_KERNEL_H_
#define SET_STRETCH_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

//{{{ Set stretch
__global__ void set_stretch_kernel(int samps, float mean, float *d_input)
{

	int t = blockIdx.x * blockDim.x + threadIdx.x;

	if (t >= 0 && t < samps)
		d_input[t] = mean;
}

//}}}

#endif

