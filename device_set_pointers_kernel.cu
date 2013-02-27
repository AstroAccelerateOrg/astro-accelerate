#ifndef POINTERS_KERNEL_H_
#define POINTERS_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include <cutil_inline.h>

//{{{ Set pointers

__global__ void set_device_pointers(int i, float **A, float *B)
{
	A[i] = B;
}

//}}}

#endif
