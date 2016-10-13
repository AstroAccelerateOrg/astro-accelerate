#ifndef SET_STRETCH_KERNEL_H_
#define SET_STRETCH_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "AstroAccelerate/params.h"

<<<<<<< HEAD
//{{{ Set stretch
=======
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
__global__ void set_stretch_kernel(int samps, float mean, float *d_input)
{

	int t  = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(t >= 0 && t <samps) d_input[t] = mean;
}

<<<<<<< HEAD
//}}}

=======
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
#endif

