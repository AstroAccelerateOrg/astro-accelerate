/* This function takes a pointer to the file pointer so that it can update the position of the file pointer
 */

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
#include "helper_cuda.h"


void allocate_memory_cpu_input(FILE **fp, size_t gpu_memory, int maxshift, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, unsigned short **input_buffer, float ****output_buffer, unsigned short **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize)
{

	*inputsize = nsamp * (size_t) nchans * sizeof(unsigned short);
	printf("\nNsamp:\t%i \tNchans:\t%i\t short:\t%i inputsize: %zu",nsamp,nchans,sizeof(unsigned short),*inputsize/1024/1024);

	*input_buffer = (unsigned short *) malloc(*inputsize);
//	cudaMallocHost((void  **) &(*input_buffer),*inputsize);
}

void allocate_memory_cpu_output(FILE **fp, size_t gpu_memory, int maxshift, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, unsigned short **input_buffer, float ****output_buffer, unsigned short **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize)
{

	*outputsize = 0;
	*output_buffer = (float ***) malloc(range * sizeof(float **));
//        checkCudaErrors(cudaMallocHost((void **) &(*output_buffer), sizeof(float **)*range));
	for (int i = 0; i < range; i++)
	{
		int total_samps = 0;
		for (int k = 0; k < num_tchunks; k++)
			total_samps += t_processed[i][k];
		//printf("\nTOTSAMPS:\t%d %d", total_samps, i);
		( *output_buffer )[i] = (float **) malloc(ndms[i] * sizeof(float *));
//	        cudaMallocHost((void ** ) &((*output_buffer)[i]), sizeof(float*)*ndms[i]);
		//if((*output_buffer)[i]) printf("\n FAILED! Could not allocate %zu bytes", ndms[i]*sizeof(float *));
		for (int j = 0; j < ndms[i]; j++)
		{
			( *output_buffer )[i][j] = (float *) malloc(( total_samps ) * sizeof(float));
//			cudaMallocHost((void **) &((*output_buffer)[i][j]),sizeof(float)*total_samps);
//			if((*output_buffer)[i][j]) printf("\n FAILED! Could not allocate %zu bytes", ndms[i]*sizeof(float *));
//			memset((*output_buffer)[i][j],0.0f,(total_samps)*sizeof(float));
		}
		*outputsize += ( total_samps ) * ndms[i] * sizeof(float);
		printf("\noutput size: %llu", (unsigned long long) sizeof( *output_buffer) );
	}
}

void allocate_memory_gpu(FILE **fp, size_t gpu_memory, int maxshift, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, unsigned short **input_buffer, unsigned short **input_buffer_small, float ****output_buffer, float **output_buffer_small, unsigned short **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize)
{

	int time_samps = NUM_STREAMS*(t_processed[0][0] + maxshift+1);
	*gpu_inputsize = (size_t) time_samps * (size_t) nchans * sizeof(unsigned short);
	printf("\n\n\n%d\n\n\n", time_samps), fflush(stdout);
	checkCudaErrors( cudaMalloc((void **) d_input, *gpu_inputsize) );
//	*input_buffer_small = (unsigned short *) malloc(*gpu_inputsize);
	checkCudaErrors( cudaMallocHost((void **) &(*input_buffer_small),*gpu_inputsize));
//		checkCudaErrors( cudaMallocManaged((void **) d_input, *gpu_inputsize) );

	if (nchans < max_ndms)
	{
		*gpu_outputsize = (size_t)time_samps * (size_t)max_ndms * sizeof(float);
	}
	else
	{
		*gpu_outputsize = (size_t)time_samps * (size_t)nchans * sizeof(float);
	}
	
	checkCudaErrors( cudaMalloc((void **) d_output, *gpu_outputsize) );
	checkCudaErrors( cudaMallocHost((void **) output_buffer_small, time_samps*sizeof(float)*max_ndms) );

	//end_t=omp_get_wtime();
	//time = (float)(end_t-start_t);
	//printf("\nGPU Malloc in: %f ", time);

	( cudaMemset(*d_output, 0, *gpu_outputsize) );

}

void allocate_memory_MSD(float **d_MSD_workarea, unsigned short **d_MSD_output_taps, float **d_MSD_interpolated, float **h_MSD_interpolated, float **h_MSD_DIT, int **gmem_peak_pos, int **temp_peak_pos, unsigned long int MSD_maxtimesamples, int MSD_DIT_widths, int nTimesamples, size_t MSD_profile_size){

//	printf("\n\n\n%lld\n\n\n", MSD_maxtimesamples), fflush(stdout);

	checkCudaErrors(cudaMalloc((void **) d_MSD_workarea, NUM_STREAMS*MSD_maxtimesamples*5.5*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **) &(*d_MSD_output_taps), NUM_STREAMS*sizeof(ushort)*2*MSD_maxtimesamples));
        checkCudaErrors(cudaMallocManaged((void **) gmem_peak_pos, NUM_STREAMS*sizeof(int)));
//        checkCudaErrors(cudaMalloc((void **) gmem_peak_pos, NUM_STREAMS*sizeof(int)));
	checkCudaErrors(cudaMallocHost((void **) temp_peak_pos, NUM_STREAMS*sizeof(int)));
	checkCudaErrors(cudaMalloc((void **) d_MSD_interpolated, NUM_STREAMS*sizeof(float)*MSD_profile_size));

	// host pinned allocation

	checkCudaErrors(cudaMallocHost((void **) h_MSD_DIT, NUM_STREAMS*180));
	checkCudaErrors(cudaMallocHost((void **) &h_MSD_interpolated, NUM_STREAMS*sizeof(float)*MSD_profile_size));

}
