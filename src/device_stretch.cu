// #include <omp.h>
#include <stdio.h>
#include "headers/params.h"
#include "device_stretch_kernel.cu"
#include "helper_cuda.h"

//{{{ Dopler Stretch 

void stretch_gpu(cudaEvent_t event, cudaStream_t stream, int acc, int samps, float tsamp, float *d_input, float *d_output)
{

	//{{{ Simple corner turn on the GPU 

	int divisions_in_t = 32;
	int num_blocks_t = samps / divisions_in_t;

	float t_zero = ( (double) tsamp ) / ( 1.0 + ( ( acc * samps * (double) tsamp ) / 599584916.0 ) );
	float multiplier = ( t_zero * acc * (double) tsamp ) / 599584916.0;
	float tsamp_inverse = 1.0 / tsamp;

//	printf("\nStretch!");
//	printf("\n%d %d", samps, acc);
//	printf("\n%d %d", divisions_in_t, num_blocks_t);

	dim3 threads_per_block(divisions_in_t);
	dim3 num_blocks(num_blocks_t);

	//double start_t, end_t;
	//start_t = omp_get_wtime();

	cudaStreamWaitEvent(stream, event, 0);
	stretch_kernel<<<num_blocks, threads_per_block, 0, stream>>>(acc, samps, tsamp, d_input, d_output, t_zero, multiplier, tsamp_inverse);
	getLastCudaError("stretch_kernel failed");
	cudaEventRecord(event, stream);

	//end_t = omp_get_wtime();
	//float time = (float)(end_t-start_t);
	//printf("\nPerformed Stretch: %f (GPU estimate)", time);
	//printf("\nCT Gops based on %.2f ops per channel per tsamp: %f",10.0,((10.0*(divisions_in_t*divisions_in_f*num_blocks_t*num_blocks_f))/(time))/1000000000.0);
	//printf("\nCT Device memory bandwidth in GB/s: %f", ((sizeof(float)+sizeof(unsigned short))*(divisions_in_t*divisions_in_f*num_blocks_t*num_blocks_f))/(time)/1000000000.0);

	//}}}
}

//}}}

