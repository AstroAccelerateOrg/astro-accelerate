#ifndef CORNERTURN_KERNEL_H_
#define CORNERTURN_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "AstroAccelerate/params.h"


//{{{ corner_turn
__global__ void simple_corner_turn_kernel(float *d_input, float *d_output, int nchans, int nsamp)
{

	int t  = blockIdx.x * blockDim.x + threadIdx.x;
	int c  = blockIdx.y * blockDim.y + threadIdx.y;
	
	d_output[c*nsamp + t] = (float)d_input[t*nchans + c];
	//if(blockIdx.x==0 && threadIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) printf("\n %d %d", nsamp, nchans);
	

}

__global__ void swap(float *d_input, float *d_output, int nchans, int nsamp)
{

	int t  = blockIdx.x * blockDim.x + threadIdx.x;
	int c  = blockIdx.y * blockDim.y + threadIdx.y;
	
	d_input[c*nsamp + t] = (float)d_output[c*nsamp + t];	

}



//}}}

#endif

