#include <cutil_inline.h>

extern "C" void load_data(float *device_pointer, float *host_pointer, size_t size, int nsamp, int maxshift);

//{{{ load_data_from_host_to_device

void load_data(float *device_pointer, float *host_pointer, size_t size, int nsamp, int maxshift) {

	//{{{ Copy data and set up the GPU constants/variables.

	//cudaEvent_t start, stop;
	//float time;
	//cudaEventCreate(&start);
	//cudaEventCreate(&stop);

	//printf("\n\tmemStart"),fflush(stdout);
	//cudaEventRecord(start,0);
	
	cutilSafeCall( cudaMemcpy(device_pointer, host_pointer, size, cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpyToSymbol("i_nsamp", &nsamp, sizeof(int)) );
	cutilSafeCall( cudaMemcpyToSymbol("i_maxshift", &maxshift, sizeof(int)) );

	//cudaEventRecord(stop, 0);
	//cudaEventSynchronize(stop);
	//cudaEventElapsedTime(&time, start, stop);
	//printf("\n\tmemStop"),fflush(stdout);
	//printf("\n\tCopied data to GPU:\t\t\t\t%lf ms", time);    
	//printf("\n\n\tEffective bandwidth in GB per second (input):\t%f\n", (((float)size)/1000000000)/(time/1000));

	//}}}
	
}

//}}}
