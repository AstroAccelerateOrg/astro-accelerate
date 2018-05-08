/* This function takes a pointer to the file pointer so that it can update the position of the file pointer
 */
 
#include <helper_cuda.h>

#include <vector_types.h>
#include <driver_functions.h>
#include <cuda_runtime.h>

// CUDA utilities and system includes
#include <vector_types.h>

#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
//#include <omp.h>
#include <cuda.h>
#include "headers/params.h"

#include "headers/device_DDTR_Plan.h"

void allocate_memory_cpu_input(unsigned short **input_buffer, size_t nsamp, size_t nchans) {

	size_t inputsize = nsamp*((size_t) nchans)*sizeof(unsigned short);
	*input_buffer = (unsigned short *) malloc(inputsize);
}

void allocate_memory_cpu_output(float ****output_buffer, DDTR_Plan *DDTR_plan) {
	int nRanges = DDTR_plan->nRanges;
	int num_tchunks = DDTR_plan->num_tchunks;
	size_t host_outputsize = 0;
	
	*output_buffer = (float ***) malloc( nRanges*sizeof(float **) );
	for (int i = 0; i < nRanges; i++) {
		size_t total_samps = 0;
		for (int k = 0; k < num_tchunks; k++)
			total_samps += DDTR_plan->t_processed[i][k];
		//printf("\nTOTSAMPS:\t%d %d", total_samps, i);
		( *output_buffer )[i] = (float **) malloc( DDTR_plan->ndms[i]*sizeof(float *));
		//if((*output_buffer)[i]) printf("\n FAILED! Could not allocate %zu bytes", ndms[i]*sizeof(float *));
		for (int j = 0; j < DDTR_plan->ndms[i]; j++) {
			( *output_buffer )[i][j] = (float *) malloc( total_samps*sizeof(float) );
			if( (*output_buffer)[i][j] == NULL ) printf("\n FAILED! Could not allocate %zu bytes", DDTR_plan->ndms[i]*sizeof(float *));
		}
		host_outputsize += total_samps*DDTR_plan->ndms[i]*sizeof(float);
		printf("\noutput size: %llu", (unsigned long long) sizeof( *output_buffer ) / 1024 / 1024 / 1024);
	}
	
	DDTR_plan->host_outputsize = host_outputsize;
}

void allocate_memory_gpu(unsigned short **d_input, float **d_output, DDTR_Plan *DDTR_plan) {
	size_t time_samps = DDTR_plan->t_processed[0][0] + (size_t) DDTR_plan->maxshift;
	printf("\n\n\n%d\n\n\n", time_samps), fflush(stdout);
	DDTR_plan->gpu_inputsize = time_samps*DDTR_plan->nchans*sizeof(unsigned short);
	checkCudaErrors( cudaMalloc((void **) d_input, DDTR_plan->gpu_inputsize) );
	printf("time_samp: %zu; input size: %zu\n", time_samps, DDTR_plan->gpu_inputsize);

	if (DDTR_plan->nchans < DDTR_plan->max_ndms) {
		DDTR_plan->gpu_outputsize = time_samps*DDTR_plan->max_ndms*sizeof(float);
	}
	else {
		DDTR_plan->gpu_outputsize = time_samps*DDTR_plan->nchans*sizeof(float);
	}
	printf("DDTR_plan->gpu_outputsize: %zu;\n", DDTR_plan->gpu_outputsize);
	checkCudaErrors( cudaMalloc((void **) d_output, DDTR_plan->gpu_outputsize) );
	checkCudaErrors( cudaMemset(*d_output, 0, DDTR_plan->gpu_outputsize) );
}



