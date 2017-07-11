//#include <omp.h>
#include <time.h>
#include <stdio.h>
#include "headers/params.h"
#include "device_corner_turn_kernel.cu"

//{{{ Corner-turn 

void corner_turn(unsigned short *d_input, float *d_output, int nchans, int nsamp)
{

	//{{{ Simple corner turn on the GPU 

	int divisions_in_t = CT;
	int divisions_in_f = CF;
	int num_blocks_t = nsamp / divisions_in_t;
	int num_blocks_f = nchans / divisions_in_f;

//	printf("\nCORNER TURN!");
//	printf("\n%d %d", nsamp, nchans);
//	printf("\n%d %d", divisions_in_t, divisions_in_f);
//	printf("\n%d %d", num_blocks_t, num_blocks_f);

	dim3 threads_per_block(divisions_in_t, divisions_in_f);
	dim3 num_blocks(num_blocks_t, num_blocks_f);

	clock_t start_t, end_t;
	start_t = clock();

	simple_corner_turn_kernel<<<num_blocks, threads_per_block>>>(d_input, d_output, nchans, nsamp);
	cudaDeviceSynchronize();
	swap<<<num_blocks, threads_per_block>>>(d_input, d_output, nchans, nsamp);
	cudaDeviceSynchronize();

	end_t = clock();
	double time = (double) ( end_t - start_t )/CLOCKS_PER_SEC;
//	printf("\nPerformed CT: %lf (GPU estimate)", time);
//	printf("\nCT Gops based on %.2f ops per channel per tsamp: %f", 10.0, ( ( 10.0 * ( divisions_in_t * divisions_in_f * num_blocks_t * num_blocks_f ) ) / ( time ) ) / 1000000000.0);
//	printf("\nCT Device memory bandwidth in GB/s: %lf", ( ( sizeof(float) + sizeof(unsigned short) ) * ( divisions_in_t * divisions_in_f * num_blocks_t * num_blocks_f ) ) / ( time ) / 1000000000.0);

	//cudaMemcpy(d_input, d_output, inputsize, cudaMemcpyDeviceToDevice);

	//}}}

}

//}}}

