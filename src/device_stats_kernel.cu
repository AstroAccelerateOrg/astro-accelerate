#ifndef STATS_KERNEL_H_
#define STATS_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "headers/params.h"

//{{{ Set stats
__global__ void stats_kernel(int half_samps, float *d_sum, float *d_sum_square, float *d_signal_power)
{

	int t = blockIdx.x * blockDim.x * STATSLOOP + threadIdx.x;

	float local = 0.0;
	float sum = 0.0;
	float sum_square = 0.0;

	for (int i = t; i < t + STATSLOOP * blockDim.x; i += blockDim.x)
	{
		local = d_signal_power[i];
		sum += local;
		sum_square += local * local;
	}
	d_sum[blockIdx.x * blockDim.x + threadIdx.x] = sum;
	d_sum_square[blockIdx.x * blockDim.x + threadIdx.x] = sum_square;
}
#endif

