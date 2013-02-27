#include <cutil_inline.h>

extern "C" void save_data(float *device_pointer, float *host_pointer, size_t size);

//{{{ save_data_from_device_to_host

void save_data(float *device_pointer, float *host_pointer, size_t size) {

	//{{{ Copy data and set up the GPU constants/variables.

//	cudaEvent_t start, stop;
//	float time;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);

//	printf("\n\tmemStart"),fflush(stdout);
//	cudaEventRecord(start,0);

	cutilSafeCall(cudaMemcpy(host_pointer, device_pointer, size,cudaMemcpyDeviceToHost));

//	cudaEventRecord(stop, 0);
//	cudaEventSynchronize(stop);
//	cudaEventElapsedTime(&time, start, stop);
//	printf("\n\tmemStop"),fflush(stdout);
//	printf("\n\tCopied data to GPU:\t\t\t\t%lf ms", time);    
//	printf("\n\n\tEffective bandwidth in GB per second (input):\t%f\n", (((float)inputsize)/1000000000)/(time/1000));

	//}}}
	
}

//}}}
