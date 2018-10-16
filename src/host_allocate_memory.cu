/* This function takes a pointer to the file pointer so that it can update the
 * position of the file pointer
 */

#include <cuda_runtime.h>
#include <driver_functions.h>
#include <vector_types.h>

// CUDA utilities and system includes
#include <vector_types.h>

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
//#include <omp.h>
#include "host_info.hpp"
#include "params.hpp"
#include <cuda.h>
#include <helper_cuda.h>

void allocate_memory_cpu_input(FILE**           fp,
                               size_t           gpu_memory,
                               size_t*          host_memory,
                               int              maxshift,
                               int              num_tchunks,
                               int              max_ndms,
                               int              total_ndms,
                               int              nsamp,
                               int              nchans,
                               int              nbits,
                               int              range,
                               int*             ndms,
                               int**            t_processed,
                               unsigned short** input_buffer,
                               float****        output_buffer,
                               unsigned short** d_input,
                               float**          d_output,
                               size_t*          gpu_inputsize,
                               size_t*          gpu_outputsize,
                               size_t*          inputsize,
                               size_t*          outputsize) {
  //	printf("\nAvailable memory: %zu MiB\n", (size_t)host_memory/1024/1024);
  *inputsize = nsamp * (size_t)nchans * sizeof(unsigned short);
  if(*host_memory < *inputsize) {
    host_mem_error((unsigned int)(*inputsize / 1024.0 / 1024.0),
                   (unsigned int)(*host_memory / 1024.0 / 1024.0),
                   "input");
  }
  *host_memory  = *host_memory - *inputsize;
  *input_buffer = (unsigned short*)malloc(*inputsize);
  //        printf("\nMemory available: %zd", *host_memory/1024/1024);
}

void allocate_memory_cpu_output(FILE**           fp,
                                size_t           gpu_memory,
                                size_t*          host_memory,
                                int              maxshift,
                                int              num_tchunks,
                                int              max_ndms,
                                int              total_ndms,
                                int              nsamp,
                                int              nchans,
                                int              nbits,
                                int              range,
                                int*             ndms,
                                int**            t_processed,
                                unsigned short** input_buffer,
                                float****        output_buffer,
                                unsigned short** d_input,
                                float**          d_output,
                                size_t*          gpu_inputsize,
                                size_t*          gpu_outputsize,
                                size_t*          inputsize,
                                size_t*          outputsize) {
  size_t estimate_outputbuffer_size = 0;
  for(int i = 0; i < range; i++)
    for(int j = 0; j < num_tchunks; j++)
      estimate_outputbuffer_size +=
          (size_t)(t_processed[i][j] * sizeof(float) * ndms[i]);
  //	size_t estimate_outputbuffer_size =
  //(size_t)(total_ndms*sizeof(float)*nsamp); 	printf("\nTotal ndms: %i; nsamp:
  //%i; mem: %zu estimated: %zu \n", total_ndms, nsamp,
  //(size_t)(*host_memory)/1024/1024, estimate_outputbuffer_size/1024/1024);
  if(*host_memory < estimate_outputbuffer_size) {
    host_mem_error((unsigned int)(estimate_outputbuffer_size / 1024.0 / 1024.0),
                   (unsigned int)(*host_memory / 1024.0 / 1024.0),
                   "output");
  }

  *outputsize    = 0;
  *output_buffer = (float***)malloc(range * sizeof(float**));
  for(int i = 0; i < range; i++) {
    int total_samps = 0;
    for(int k = 0; k < num_tchunks; k++)
      total_samps += t_processed[i][k];
    // printf("\nTOTSAMPS:\t%d %d", total_samps, i);
    (*output_buffer)[i] = (float**)malloc(ndms[i] * sizeof(float*));
    // if((*output_buffer)[i]) printf("\n FAILED! Could not allocate %zu bytes",
    // ndms[i]*sizeof(float *));
    for(int j = 0; j < ndms[i]; j++) {
      (*output_buffer)[i][j] = (float*)malloc((total_samps) * sizeof(float));
      // if((*output_buffer)[i][j]) printf("\n FAILED! Could not allocate %zu
      // bytes", ndms[i]*sizeof(float *));
      //			memset((*output_buffer)[i][j],0.0f,(total_samps)*sizeof(float));
    }
    *outputsize += (total_samps)*ndms[i] * sizeof(float);
  }
  *host_memory = *host_memory - *outputsize;
  //      printf("\noutput: Memory available: %zu", (*host_memory)/1024/1024);
  //	printf("\noutput size: %llu",(unsigned long
  //long)(*outputsize/1024/1024));
}

void allocate_memory_gpu(FILE**           fp,
                         size_t           gpu_memory,
                         int              maxshift,
                         int              num_tchunks,
                         int              max_ndms,
                         int              total_ndms,
                         int              nsamp,
                         int              nchans,
                         int              nbits,
                         int              range,
                         int*             ndms,
                         int**            t_processed,
                         unsigned short** input_buffer,
                         float****        output_buffer,
                         unsigned short** d_input,
                         float**          d_output,
                         size_t*          gpu_inputsize,
                         size_t*          gpu_outputsize,
                         size_t*          inputsize,
                         size_t*          outputsize) {

  int time_samps = t_processed[0][0] + maxshift;
  printf("\n\n\n%d\n\n\n", time_samps), fflush(stdout);
  *gpu_inputsize = (size_t)time_samps * (size_t)nchans * sizeof(unsigned short);
  checkCudaErrors(cudaMalloc((void**)d_input, *gpu_inputsize));

  if(nchans < max_ndms) {
    *gpu_outputsize = (size_t)time_samps * (size_t)max_ndms * sizeof(float);
  } else {
    *gpu_outputsize = (size_t)time_samps * (size_t)nchans * sizeof(float);
  }
  checkCudaErrors(cudaMalloc((void**)d_output, *gpu_outputsize));

  // end_t=omp_get_wtime();
  // time = (float)(end_t-start_t);
  // printf("\nGPU Malloc in: %f ", time);

  (cudaMemset(*d_output, 0, *gpu_outputsize));
}

void allocate_memory_MSD(float**           d_MSD_workarea,
                         unsigned short**  d_MSD_output_taps,
                         float**           d_MSD_interpolated,
                         unsigned long int MSD_maxtimesamples,
                         int               MSD_DIT_widths,
                         int               nTimesamples,
                         size_t            MSD_profile_size) {

  checkCudaErrors(cudaMalloc((void**)d_MSD_workarea,
                             MSD_maxtimesamples * 5.5 * sizeof(float)));
  checkCudaErrors(cudaMalloc((void**)&(*d_MSD_output_taps),
                             sizeof(ushort) * 2 * MSD_maxtimesamples));
  //        checkCudaErrors(cudaMallocManaged((void **) gmem_peak_pos,
  //        sizeof(int)));
  //       checkCudaErrors(cudaMallocHost((void **) temp_peak_pos,
  //       sizeof(int)));
  checkCudaErrors(
      cudaMalloc((void**)d_MSD_interpolated, sizeof(float) * MSD_profile_size));
}
