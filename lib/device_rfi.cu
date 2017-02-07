
#include <omp.h>
#include <stdio.h>
#include "AstroAccelerate/params.h"
#include "device_rfi_kernel.cu"

//{{{ rfi_gpu

void rfi_gpu(unsigned short *d_input, int nchans, int nsamp) {

	int divisions_in_f  = 32;
	int num_blocks_f    = nchans/divisions_in_f;

	printf("\nCORNER TURN!");
	printf("\n%d %d", nsamp, nchans);
	printf("\n%d %d", divisions_in_f, 1);
	printf("\n%d %d", num_blocks_f, 1);

	dim3 threads_per_block(divisions_in_f, 1);
	dim3 num_blocks(num_blocks_f,1);

	double start_t, end_t;
	start_t = omp_get_wtime();

	rfi_gpu_kernel<<< num_blocks, threads_per_block >>>(d_input, nchans, nsamp);
	cudaDeviceSynchronize();

	end_t = omp_get_wtime();
	float time = (float)(end_t-start_t);
	printf("\nPerformed RFI: %f (GPU estimate)", time);

	//}}}

}

//}}}

