#ifndef CORNERTURN_KERNEL_H_
#define CORNERTURN_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "AstroAccelerate/params.h"

<<<<<<< HEAD

//{{{ corner_turn
__global__ void simple_corner_turn_kernel(unsigned short *d_input, float *d_output, int nchans, int nsamp)
{

=======
__global__ void simple_corner_turn_kernel(unsigned short *d_input, float *d_output, int nchans, int nsamp)
{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	int t  = blockIdx.x * blockDim.x + threadIdx.x;
	int c  = blockIdx.y * blockDim.y + threadIdx.y;
	
	d_output[c*nsamp + t] = (float)__ldg(&d_input[t*nchans + c]);
<<<<<<< HEAD

=======
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
}

__global__ void swap(unsigned short *d_input, float *d_output, int nchans, int nsamp)
{
<<<<<<< HEAD

=======
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	int t  = blockIdx.x * blockDim.x + threadIdx.x;
	int c  = blockIdx.y * blockDim.y + threadIdx.y;
	
	d_input[c*nsamp + t] = (unsigned short)__ldg(&d_output[c*nsamp + t]);	
<<<<<<< HEAD

}



//}}}

#endif

=======
}

#endif
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
