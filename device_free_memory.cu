#include <cutil_inline.h>

extern "C" void free_device_memory(float *device_pointer);

//{{{ init_gpu

void free_device_memory(float *device_pointer) {

	//{{{ Free the memory

//	cudaEvent_t start, stop;
//	float time;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);

//	printf("\n\n\tfreeStart"),fflush(stdout);	
//	cudaEventRecord(start,0);
	
	cutilSafeCall( cudaFree(device_pointer));

//	cudaEventRecord(stop, 0);
//	cudaEventSynchronize(stop);
//	cudaEventElapsedTime(&time, start, stop);
//	printf("\n\tfreeStart"),fflush(stdout);	
//	printf("\n\tGPU memory free:\t\t\t\t%lf ms", time);

//	cudaEventDestroy(start); 
//	cudaEventDestroy(stop);

	//}}}
	
}

//}}}

