
//#include <omp.h>
#include <time.h>
#include <stdio.h>
#include "headers/params.h"
#include "device_zero_dm_outliers_kernel.cu"

//{{{ zero_dm

void zero_dm_outliers(unsigned short *d_input, int nchans, int nsamp) {

	int divisions_in_t  = CT;
	int num_blocks_t    = nsamp/divisions_in_t;

	printf("\nCORNER TURN!");
	printf("\n%d %d", nsamp, nchans);
	printf("\n%d %d", divisions_in_t, 1);
	printf("\n%d %d", num_blocks_t, 1);

	dim3 threads_per_block(divisions_in_t, 1);
	dim3 num_blocks(num_blocks_t,1);

	clock_t start_t, end_t;
	start_t = clock();

	zero_dm_outliers_kernel<<< num_blocks, threads_per_block >>>(d_input, nchans, nsamp);
	cudaDeviceSynchronize();


	end_t = clock();
	double time = (double)(end_t-start_t) / CLOCKS_PER_SEC;
	printf("\nPerformed ZDM: %lf (GPU estimate)", time);

	//}}}

}

//}}}

