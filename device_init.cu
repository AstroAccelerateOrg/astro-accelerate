#include <cutil_inline.h>

extern "C" void init_gpu(int nchans, float *dmshifts);

//{{{ init_gpu

void init_gpu(int nchans, float *dmshifts) {

	int i, devCount, device;

	//{{{ Determine what GPUs are availible
	
	cudaGetDeviceCount(&devCount);
	printf("\n\nCUDA Device Query...");
	printf("\nThere are %d CUDA devices.", devCount);
	/*
	for(i = 0; i < devCount; ++i) {
		
		// Get device properties
		printf("\n\nCUDA Device #%d\n", i);

		cudaDeviceProp devProp;
	        cudaGetDeviceProperties(&devProp, i);
		printf("\nMajor revision number:         %d",  devProp.major);
		printf("\nMinor revision number:         %d",  devProp.minor);
		printf("\nName:                          %s",  devProp.name);
		printf("\nTotal global memory:           %u",  devProp.totalGlobalMem);
		printf("\nTotal shared memory per block: %u",  devProp.sharedMemPerBlock);
		printf("\nTotal registers per block:     %d",  devProp.regsPerBlock);
		printf("\nWarp size:                     %d",  devProp.warpSize);
		printf("\nMaximum memory pitch:          %u",  devProp.memPitch);
		printf("\nMaximum threads per block:     %d",  devProp.maxThreadsPerBlock);
		for (int i = 0; i < 3; ++i) printf("\nMaximum dimension %d of block:  %d", i, devProp.maxThreadsDim[i]);
		for (int i = 0; i < 3; ++i) printf("\nMaximum dimension %d of grid:   %d", i, devProp.maxGridSize[i]);
		printf("\nClock rate:                    %d",  devProp.clockRate);
		printf("\nTotal constant memory:         %u",  devProp.totalConstMem);
		printf("\nTexture alignment:             %u",  devProp.textureAlignment);
		printf("\nConcurrent copy and execution: %s",  (devProp.deviceOverlap ? "Yes" : "No"));
		printf("\nNumber of multiprocessors:     %d",  devProp.multiProcessorCount);
		printf("\nKernel execution timeout:      %s",  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
	}
	*/

	//}}}

	//{{{ Pick a Tesla GPU card and set the device flags

	printf("\n\n\tinitStart"),fflush(stdout);

	for(i = 0; i < devCount; ++i) {

		cudaDeviceProp devProp;
	        cudaGetDeviceProperties(&devProp, i);

		if(!strncmp("Tesla", devProp.name, 5)) device = i;
	}

	printf("\n\tUsing device:\t\t\t%d", device);
	
	cudaEvent_t start, stop;
	float time;
	
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start,0);

	cudaSetDevice(device);
	cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	printf("\n\tinitStop"),fflush(stdout);
	printf("\n\tInitialised GPU:\t\t%lf ms\n", time);    

	//}}}

	//{{{ Copy data and set up the GPU constants/variables.
	
	printf("\n\tnchans:\t\t\t\t%d", nchans);
	printf("\n\tmemStart"),fflush(stdout);
	cudaEventRecord(start,0);
	
	cutilSafeCall( cudaMemcpyToSymbol(dm_shifts, dmshifts,nchans * sizeof(float)) );
	cutilSafeCall( cudaMemcpyToSymbol("i_nchans", &nchans, sizeof(int)) );


	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	printf("\n\tmemStop"),fflush(stdout);
	printf("\n\tCopied data to GPU:\t\t%lf ms", time);    

	//}}}
	
}

//}}}

