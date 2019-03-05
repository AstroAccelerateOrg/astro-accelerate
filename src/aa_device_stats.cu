#include <time.h>
#include <stdio.h>
#include <cufft.h>
#include "aa_params.hpp"
#include "aa_device_stats_kernel.hpp"
#include "aa_log.hpp"


namespace astroaccelerate {
  /** \brief Returns stats. */
  void stats_gpu(cudaEvent_t event, cudaStream_t stream, int samps, float *mean, float *stddev, float *h_signal_power, float *d_signal_power) {
    int a;
    int trials = ( 2 * ACCMAX + ACCSTEP ) / ACCSTEP;
    
    int half_samps = samps / 2;
    int acc_size = half_samps * trials;

    int divisions = STATST;
    int blocks = (int) floor((float) acc_size / divisions / STATSLOOP);

    dim3 threads_per_block(divisions);
    dim3 num_blocks(blocks);

    int size = (int) floor((float) acc_size / STATSLOOP);

    float* d_sum;
    cudaError_t e = cudaMalloc((void** )&d_sum, size * sizeof(float));

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMalloc in aa_device_stats.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    float* d_sum_square;
    e = cudaMalloc((void** )&d_sum_square, size * sizeof(float));

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMalloc in aa_device_stats.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    float* h_sum;
    e = cudaMallocHost((void** )&h_sum, size * sizeof(float));

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMallocHost in aa_device_stats.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    float* h_sum_square;
    e = cudaMallocHost((void** )&h_sum_square, size * sizeof(float));

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMallocHost in aa_device_stats.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    cudaStreamWaitEvent(stream, event, 0);
    call_kernel_stats_kernel(num_blocks, threads_per_block, 0, stream, half_samps, d_sum, d_sum_square, d_signal_power);
    //getLastCudaError("power_kernel failed");
    cudaEventRecord(event, stream);

    cudaStreamWaitEvent(stream, event, 0);
    e = cudaMemcpyAsync(h_sum, d_sum, size * sizeof(float), cudaMemcpyDeviceToHost, stream);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpyAsync in aa_device_stats.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    
    e = cudaMemcpyAsync(h_sum_square, d_sum_square, size * sizeof(float), cudaMemcpyDeviceToHost, stream);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpyAsync in aa_device_stats.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    cudaEventRecord(event, stream);
    cudaStreamSynchronize(stream);

    float total_sum = 0.0;
    float total_sum_square = 0.0;

    for (a = 0; a < size; a++)
      {
	total_sum += ( h_sum[a] );
	total_sum_square += ( h_sum_square[a] );
      }
    *mean = (float) ( total_sum / ( acc_size ) );  // Mean for data sample
    *stddev = (float) sqrt(( total_sum_square - acc_size * ( *mean ) * ( *mean ) ) / ( acc_size - 1 ));

    cudaFree(d_sum);
    cudaFree(d_sum_square);
    cudaFreeHost(h_sum);
    cudaFreeHost(h_sum_square);

  }
} //namespace astroaccelerate
