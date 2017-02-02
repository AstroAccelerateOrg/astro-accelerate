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
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, dev);

	printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

	// Console log
	cudaDriverGetVersion(&driverVersion);
	cudaRuntimeGetVersion(&runtimeVersion);
	printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n", driverVersion / 1000, ( driverVersion % 100 ) / 10, runtimeVersion / 1000, ( runtimeVersion % 100 ) / 10);
	printf("  CUDA Capability Major/Minor version number:    %d.%d\n", deviceProp.major, deviceProp.minor);
	printf("  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n", (float) deviceProp.totalGlobalMem / 1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);
	printf("  GPU Clock rate:                                %.0f MHz (%0.2f GHz)\n", deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);
	printf("  Memory Clock rate:                             %.0f Mhz\n", deviceProp.memoryClockRate * 1e-3f);
	printf("  Memory Bus Width:                              %d-bit\n", deviceProp.memoryBusWidth);
	if (deviceProp.l2CacheSize)
	{
		printf("  L2 Cache Size:                                 %d bytes\n", deviceProp.l2CacheSize);
	}
	printf("  Maximum Texture Dimension Size (x,y,z)         1D=(%d), 2D=(%d, %d), 3D=(%d, %d, %d)\n", deviceProp.maxTexture1D, deviceProp.maxTexture2D[0], deviceProp.maxTexture2D[1], deviceProp.maxTexture3D[0], deviceProp.maxTexture3D[1], deviceProp.maxTexture3D[2]);
	printf("  Maximum Layered 1D Texture Size, (num) layers  1D=(%d), %d layers\n", deviceProp.maxTexture1DLayered[0], deviceProp.maxTexture1DLayered[1]);
	printf("  Maximum Layered 2D Texture Size, (num) layers  2D=(%d, %d), %d layers\n", deviceProp.maxTexture2DLayered[0], deviceProp.maxTexture2DLayered[1], deviceProp.maxTexture2DLayered[2]);
	printf("  Total amount of constant memory:               %lu bytes\n", deviceProp.totalConstMem);
	printf("  Total amount of shared memory per block:       %lu bytes\n", deviceProp.sharedMemPerBlock);
	printf("  Total number of registers available per block: %d\n", deviceProp.regsPerBlock);
	printf("  Warp size:                                     %d\n", deviceProp.warpSize);
	printf("  Maximum number of threads per multiprocessor:  %d\n", deviceProp.maxThreadsPerMultiProcessor);
	printf("  Maximum number of threads per block:           %d\n", deviceProp.maxThreadsPerBlock);
	printf("  Max dimension size of a thread block (x,y,z): (%d, %d, %d)\n", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
	printf("  Max dimension size of a grid size    (x,y,z): (%d, %d, %d)\n", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
	printf("  Maximum memory pitch:                          %lu bytes\n", deviceProp.memPitch);
	printf("  Texture alignment:                             %lu bytes\n", deviceProp.textureAlignment);
	printf("  Concurrent copy and kernel execution:          %s with %d copy engine(s)\n", (
	    deviceProp.deviceOverlap ? "Yes" : "No" ), deviceProp.asyncEngineCount);
	printf("  Run time limit on kernels:                     %s\n",
	    deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
	printf("  Integrated GPU sharing Host Memory:            %s\n",
	    deviceProp.integrated ? "Yes" : "No");
	printf("  Support host page-locked memory mapping:       %s\n",
	    deviceProp.canMapHostMemory ? "Yes" : "No");
	printf("  Alignment requirement for Surfaces:            %s\n",
	    deviceProp.surfaceAlignment ? "Yes" : "No");
	printf("  Device has ECC support:                        %s\n",
	    deviceProp.ECCEnabled ? "Enabled" : "Disabled");
	printf("  CUDA Device Driver Mode (TCC or WDDM):         %s\n",
	    deviceProp.tccDriver ? "TCC (Tesla Compute Cluster Driver)" : "WDDM (Windows Display Driver Model)");
	printf("  Device supports Unified Addressing (UVA):      %s\n",
	    deviceProp.unifiedAddressing ? "Yes" : "No");
	printf("  Device PCI Bus ID / PCI location ID:           %d / %d\n", deviceProp.pciBusID, deviceProp.pciDeviceID);

	size_t free, total;

	cudaMemGetInfo(&free, &total);
	*gpu_memory = ( free );
}

//}}}

