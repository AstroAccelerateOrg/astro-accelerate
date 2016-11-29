#include <stdio.h>
#include "AstroAccelerate/params.h"

// CUDA-C includes
#include <cuda.h>
#include <cuda_runtime.h>

//extern "C" void init_gpu(int argc, char **argv, int enable_debug, size_t *gpu_memory);

//{{{ init_gpu

void init_gpu(int argc, char **arg, int enable_debug, size_t *gpu_memory)
{

	int deviceCount = 0;
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

	if (error_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int) error_id, cudaGetErrorString(error_id));
		printf("Result = FAIL\n");
		exit(EXIT_FAILURE);
	}

	// This function call returns 0 if there are no CUDA capable devices.
	if (deviceCount == 0)
	{
		printf("There are no available device(s) that support CUDA\n");
	}
	else
	{
		printf("Detected %d CUDA Capable device(s)\n", deviceCount);
	}

	int dev, driverVersion = 0, runtimeVersion = 0;
	dev = CARD;

	cudaSetDevice(dev);
	size_t free, total;

	cudaMemGetInfo(&free, &total);
	*gpu_memory = ( free/4 );
}

//}}}

