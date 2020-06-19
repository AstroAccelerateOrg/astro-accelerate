#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include "aa_log.hpp"

namespace astroaccelerate {
  /** \brief Copy data and set up the GPU constants/variables. */
  void save_data(float *device_pointer, float *host_pointer, size_t size) {
    cudaMemcpy(host_pointer, device_pointer, size, cudaMemcpyDeviceToHost);
  }

  /** \brief Copy data and set up the GPU constants/variables. */
  void save_data_offset(float *device_pointer, size_t device_offset, float *host_pointer, size_t host_offset, size_t size) {
	cudaMemcpy(host_pointer + host_offset, device_pointer + device_offset, size, cudaMemcpyDeviceToHost);
  }

  /** \brief Copy data and set up the GPU constants/variables. */
  void save_data_offset_stream(int dm_range, int current_time_chunk, int **t_processed, long int inc, int *inBin, int const*const ndms,
		  		float *d_DDTR_output, float ***output_buffer){
		  
	// the number of threads and streams needs to be same value
	int nStreams = 16;
	cudaStream_t stream_copy[16];
	cudaError_t e;

	float *h_array_pinned;
	size_t data_size = (size_t)(sizeof(float) * (size_t) t_processed[dm_range][current_time_chunk]);
	e = cudaMallocHost((void **) &h_array_pinned, data_size*nStreams);
	if (e != cudaSuccess) {
		LOG(log_level::error, "Could not create pinned memory on host for DDTR D2H copy (" + std::string(cudaGetErrorString(e)) + ")");
	}

	for (int i = 0; i < nStreams; i++){
        	e = cudaStreamCreate(&stream_copy[i]);
		if (e != cudaSuccess) {
			LOG(log_level::error, "Could not create stream in DDTR D2H copy (" + std::string(cudaGetErrorString(e)) + ")");
		}
	}

	#pragma omp parallel for num_threads(nStreams) shared(h_array_pinned, output_buffer, d_DDTR_output, dm_range, data_size, stream_copy)
	for (size_t k = 0; k < (size_t) ndms[dm_range]; k++) {
		int id_stream = omp_get_thread_num();
		size_t device_offset = (size_t) (k * (size_t) t_processed[dm_range][current_time_chunk]);
		size_t host_offset = (size_t) (inc / inBin[dm_range]);
		size_t stream_offset = (size_t) (t_processed[dm_range][current_time_chunk]);
		cudaMemcpyAsync(h_array_pinned + id_stream*stream_offset, d_DDTR_output + device_offset, data_size, cudaMemcpyDeviceToHost, stream_copy[id_stream]);
		cudaStreamSynchronize(stream_copy[id_stream]);
		memcpy(output_buffer[dm_range][k] + host_offset, h_array_pinned + id_stream*stream_offset, data_size);
	}

	for (int i = 0; i < nStreams; i++){
        	e = cudaStreamDestroy(stream_copy[i]);
		if (e != cudaSuccess) {
			LOG(log_level::error, "Could not destroy stream in DDTR D2H copy (" + std::string(cudaGetErrorString(e)) + ")");
		}
	}
	e = cudaFreeHost(h_array_pinned);
	if (e != cudaSuccess){
		LOG(log_level::error, "Could not free host pinned memory in DDTR D2H copy (" + std::string(cudaGetErrorString(e)) + ")");
	}
  }
  

} //namespace astroaccelerate
