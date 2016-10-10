#include <omp.h>
#include <stdio.h>
#include "AstroAccelerate/params.h"
#include "device_corner_turn_kernel.cu"

void corner_turn(unsigned short *d_input, float *d_output, int nchans, int nsamp)
{

	int divisions_in_t = CT;
	int divisions_in_f = CF;
	int num_blocks_t   = nsamp/divisions_in_t;
	int num_blocks_f   = nchans/divisions_in_f;

	printf("\nCORNER TURN!");
	printf("\n%d %d", nsamp, nchans);
	printf("\n%d %d", divisions_in_t, divisions_in_f);
	printf("\n%d %d", num_blocks_t, num_blocks_f);

	dim3 threads_per_block(divisions_in_t, divisions_in_f);
	dim3 num_blocks(num_blocks_t,num_blocks_f);

	double start_t, end_t;
	start_t = omp_get_wtime();

	simple_corner_turn_kernel<<< num_blocks, threads_per_block >>>(d_input, d_output, nchans, nsamp);
	cudaDeviceSynchronize();
	swap<<< num_blocks, threads_per_block >>>(d_input, d_output, nchans, nsamp);
	cudaDeviceSynchronize();

	end_t = omp_get_wtime();
	float time = (float)(end_t-start_t);
	printf("\nPerformed CT: %f (GPU estimate)", time);
	printf("\nCT Gops based on %.2f ops per channel per tsamp: %f",10.0,((10.0*(divisions_in_t*divisions_in_f*num_blocks_t*num_blocks_f))/(time))/1000000000.0);
	printf("\nCT Device memory bandwidth in GB/s: %f", ((sizeof(float)+sizeof(unsigned short))*(divisions_in_t*divisions_in_f*num_blocks_t*num_blocks_f))/(time)/1000000000.0);

}