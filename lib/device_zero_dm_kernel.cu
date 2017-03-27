#ifndef ZERODM_KERNEL_H_
#define ZERODM_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"


//{{{ zero dm kernel - needs cleaning and optimizing // WA 21/10/16
__global__ void zero_dm_kernel(unsigned short *d_input, int nchans, int nsamp)
{

	int t  = blockIdx.x * blockDim.x + threadIdx.x;

	float sum = 0.0f;
	for(int c = 0; c < nchans; c++) sum+=(float)__ldg(&d_input[t*nchans + c]);
	sum = (sum/(float)nchans);
	for(int c = 0; c < nchans; c++) d_input[t*nchans + c]=(unsigned short)((unsigned char)((float)d_input[t*nchans + c]-sum));
}

//}}}

#endif

