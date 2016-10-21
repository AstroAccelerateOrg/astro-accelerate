
#include <omp.h>
#include <stdio.h>
#include "AstroAccelerate/params.h"
#include "device_zero_dm_kernel.cu"

//{{{ zero_dm

void zero_dm(unsigned short *d_input, int nchans, int nsamp) {

	int divisions_in_t  = CT;
	int num_blocks_t    = nsamp/divisions_in_t;

	printf("\nCORNER TURN!");
	printf("\n%d %d", nsamp, nchans);
	printf("\n%d %d", divisions_in_t, 1);
	printf("\n%d %d", num_blocks_t, 1);

	dim3 threads_per_block(divisions_in_t, 1);
	dim3 num_blocks(num_blocks_t,1);

	double start_t, end_t;
	start_t = omp_get_wtime();

	zero_dm_kernel<<< num_blocks, threads_per_block >>>(d_input, nchans, nsamp);
	cudaDeviceSynchronize();

	end_t = omp_get_wtime();
	float time = (float)(end_t-start_t);
	printf("\nPerformed ZDM: %f (GPU estimate)", time);

	//}}}

}

//}}}

