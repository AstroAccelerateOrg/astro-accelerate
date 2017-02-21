//#include <omp.h>
#include <stdio.h>
#include <cufft.h>
#include "AstroAccelerate/params.h"
#include "device_power_kernel.cu"
#include "helper_cuda.h"

//{{{ Dopler Stretch 

void power_gpu(cudaEvent_t event, cudaStream_t stream, int samps, int acc, cufftComplex *d_signal_fft, float *d_signal_power)
{

	int half_samps = samps / 2;

	int divisions_in_t = 32;
	int num_blocks_t = half_samps / divisions_in_t;

	dim3 threads_per_block(divisions_in_t);
	dim3 num_blocks(num_blocks_t);

	cudaStreamWaitEvent(stream, event, 0);
	power_kernel<<<num_blocks, threads_per_block, 0, stream>>>(half_samps, acc, d_signal_fft, d_signal_power);
	getLastCudaError("power_kernel failed");
	cudaEventRecord(event, stream);
}

//}}}

