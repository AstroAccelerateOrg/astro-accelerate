#include <cutil_inline.h>

extern "C" float *malloc_gpu(size_t size, int zero_mem);

//{{{ malloc_gpu

float *malloc_gpu(size_t size, int zero_mem) {

	//{{{ Allocate GPU memory
	
//	cudaEvent_t start, stop;
//	float time;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);

//	printf("\n\n\tallocStart"),fflush(stdout);
//	cudaEventRecord(start,0);

	float *device_pointer;

	cutilSafeCall( cudaMalloc((void **) &device_pointer, size));

	if(zero_mem == 0) cutilSafeCall( cudaMemset(device_pointer, 0, size));

//	cudaEventRecord(stop, 0);
//	cudaEventSynchronize(stop);
//	cudaEventElapsedTime(&time, start, stop);
//	printf("\n\tallocStop"),fflush(stdout);        
//	printf("\n\tAllocated memory on GPU:\t\t\t%lf ms\n", time);

	//}}}
	
	return device_pointer;
	
}

//}}}
