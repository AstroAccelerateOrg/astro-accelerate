#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <time.h>
#include <string>

#include "aa_log.hpp"
#include "aa_params.hpp"
#include "aa_device_stats.hpp"
#include "aa_device_stretch.hpp"
#include "aa_device_set_stretch.hpp"
#include "aa_device_power.hpp"

namespace astroaccelerate {

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
  inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
    if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
    }
  }

  /**
   * \brief Example FFT.
   * \todo Clarify the difference between this function and other fdas functions.
   */
  void acceleration(int range, int nsamp, int max_ndms, int processed, int nboots, int num_trial_bins, int navdms, float narrow, float wide, int nsearch, float aggression, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp) {
    printf("\n");

    printf("[1DCUFFT] is starting...\n");

    size_t size;
    int a; //j;
    float mean, stddev;

    for (int i = 0; i < range; i++) {

      //double total = 0.0;

      cudaStream_t stream_e;
      //cudaError_t result_e;
      gpuErrchk(cudaStreamCreate(&stream_e));

      cudaEvent_t event_e;
      cudaEventCreate(&event_e);

      cudaStream_t stream_o;
      //cudaError_t result_o;
      gpuErrchk(cudaStreamCreate(&stream_o));

      cudaEvent_t event_o;
      cudaEventCreate(&event_o);

      int samps = processed / inBin[i];

      printf("\nsamps:\t%d", samps);
      int nearest = (int) floorf(log2f((float) samps));
      printf("\nnearest:\t%d", nearest);
      samps = (int) powf(2.0, nearest);
      printf("\nsamps:\t%d", samps);

      // Allocate memory for signal even
      float* d_signal_in_e;
      size = samps * sizeof(float);
      printf("\nSize of GPU input signal:\t%zu MB", size / 1024 / 1024);
      cudaError_t e = cudaMalloc((void** )&d_signal_in_e, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float* d_signal_transformed_e;
      size = samps * sizeof(float);
      printf("\nSize of GPU stretched signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMalloc((void** )&d_signal_transformed_e, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      cufftComplex* d_signal_fft_e;
      size = ( samps / 2 + 1 ) * sizeof(cufftComplex);
      printf("\nSize of GPU output signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMalloc((void** )&d_signal_fft_e, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float* d_signal_power_e;
      size = sizeof(float) * ( samps / 2 ) * ( 2 * ACCMAX + ACCSTEP ) / ACCSTEP;
      printf("\nSize of GPU power signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMalloc((void** )&d_signal_power_e, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float2* h_signal_e;
      size = ( samps ) * sizeof(float2);
      printf("\nSize of host output signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMallocHost((void** )&h_signal_e, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float* h_signal_transformed_e;
      size = samps * sizeof(float);
      printf("\nSize of GPU stretched signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMallocHost((void** )&h_signal_transformed_e, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float* h_signal_power_e;
      size = sizeof(float) * ( samps / 2 ) * ( 2 * ACCMAX + ACCSTEP ) / ACCSTEP;
      printf("\nSize of total host power signal:\t%zu MB", size / 1024 / 1024), fflush(stdout);
      e = cudaMallocHost((void** )&h_signal_power_e, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      // Allocate memory for signal odd
      float* d_signal_in_o;
      size = samps * sizeof(float);
      printf("\nSize of GPU input signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMalloc((void** )&d_signal_in_o, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float* d_signal_transformed_o;
      size = samps * sizeof(float);
      printf("\nSize of GPU stretched signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMalloc((void** )&d_signal_transformed_o, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      cufftComplex* d_signal_fft_o;
      size = ( samps / 2 + 1 ) * sizeof(cufftComplex);
      printf("\nSize of GPU output signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMalloc((void** )&d_signal_fft_o, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float* d_signal_power_o;
      size = sizeof(float) * ( samps / 2 ) * ( 2 * ACCMAX + ACCSTEP ) / ACCSTEP;
      printf("\nSize of GPU power signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMalloc((void** )&d_signal_power_o, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float2* h_signal_o;
      size = ( samps ) * sizeof(float2);
      printf("\nSize of host output signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMallocHost((void** )&h_signal_o, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float* h_signal_transformed_o;
      size = samps * sizeof(float);
      printf("\nSize of GPU stretched signal:\t%zu MB", size / 1024 / 1024);
      e = cudaMallocHost((void** )&h_signal_transformed_o, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      float* h_signal_power_o;
      size = sizeof(float) * ( samps / 2 ) * ( 2 * ACCMAX + ACCSTEP ) / ACCSTEP;
      printf("\nSize of total host power signal:\t%zu MB", size / 1024 / 1024), fflush(stdout);
      e = cudaMallocHost((void** )&h_signal_power_o, size);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }

      // CUFFT plan even
      cufftHandle plan_e;
      cufftPlan1d(&plan_e, samps, CUFFT_R2C, 1);
      cufftSetStream(plan_e, stream_e);

      // CUFFT plan odd
      cufftHandle plan_o;
      cufftPlan1d(&plan_o, samps, CUFFT_R2C, 1);
      cufftSetStream(plan_o, stream_o);

      int trials = ( 2 * ACCMAX + ACCSTEP ) / ACCSTEP;

      // Transfer even memory asynchronously
      //TEST:checkCudaErrors(cudaMemcpyAsync(d_signal_in_e, output_buffer[i][230],   samps*sizeof(float), cudaMemcpyHostToDevice, stream_e));
      e = cudaMemcpyAsync(d_signal_in_e, output_buffer[i][0], samps * sizeof(float), cudaMemcpyHostToDevice, stream_e);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      cudaEventRecord(event_e, stream_e);

      // Cacluclate even dm
      for (a = 0; a < trials; a++)
	{
	  int acc = -ACCMAX + a * ACCSTEP;
	  float mean = 127.959f;
	  set_stretch_gpu(event_e, stream_e, samps, mean, d_signal_transformed_e);
	  stretch_gpu(event_e, stream_e, acc, samps, tsamp, d_signal_in_e, d_signal_transformed_e);
	  cudaStreamWaitEvent(stream_e, event_e, 0);
	  cufftResult e = cufftExecR2C(plan_e, (float * )d_signal_transformed_e, (cufftComplex * )d_signal_fft_e);

	  if(e != CUFFT_SUCCESS) {
	    LOG(log_level::error, "Could not cufftExecR2C in aa_host_acceleration.cu");
	  }
	  
	  power_gpu(event_e, stream_e, samps, a, d_signal_fft_e, d_signal_power_e);
	}

      for (int dm_count = 1; dm_count < ndms[i] - 1; dm_count += 2)
	{
	  cudaStreamWaitEvent(stream_o, event_o, 0);
	  cudaError_t e = cudaMemcpyAsync(d_signal_in_o, output_buffer[i][dm_count], samps * sizeof(float), cudaMemcpyHostToDevice, stream_o);
	  
	  if(e != cudaSuccess) {
	    LOG(log_level::error, "Could not cudaMemcpyAsync in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
	  }
	  
	  cudaEventRecord(event_o, stream_o);

	  // Cacluclate odd dm
	  for (a = 0; a < trials; a++)
	    {
	      int acc = -ACCMAX + a * ACCSTEP;
	      float mean = 127.959f;
	      set_stretch_gpu(event_o, stream_o, samps, mean, d_signal_transformed_o);
	      stretch_gpu(event_o, stream_o, acc, samps, tsamp, d_signal_in_o, d_signal_transformed_o);
	      cufftResult cufft_e = cufftExecR2C(plan_o, (float * )d_signal_transformed_o, (cufftComplex * )d_signal_fft_o);
	      
	      if(cufft_e != CUFFT_SUCCESS) {
		LOG(log_level::error, "Could not cufftExecR2C in aa_host_acceleration.cu");
	      }
	      
	      cudaStreamWaitEvent(stream_o, event_o, 0);
	      power_gpu(event_o, stream_o, samps, a, d_signal_fft_o, d_signal_power_o);
	    }

	  // Threshold even f-fdot plane
	  cudaStreamSynchronize(stream_e);
	  stats_gpu(event_e, stream_e, samps, &mean, &stddev, h_signal_power_e, d_signal_power_e);
	    
	  e = cudaMemcpyAsync(d_signal_in_e, output_buffer[i][dm_count + 1], samps * sizeof(float), cudaMemcpyHostToDevice, stream_e);
	  
	  if(e != cudaSuccess) {
	    LOG(log_level::error, "Could not cudaMemcpyAsync in aa_host_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
	  }
	  
	  cudaEventRecord(event_e, stream_e);

	  // Cacluclate even dm
	  for (a = 0; a < trials; a++)
	    {
	      int acc = -ACCMAX + a * ACCSTEP;
	      float mean = 127.959f;
	      set_stretch_gpu(event_e, stream_e, samps, mean, d_signal_transformed_e);
	      stretch_gpu(event_e, stream_e, acc, samps, tsamp, d_signal_in_e, d_signal_transformed_e);
	      cudaStreamWaitEvent(stream_e, event_e, 0);
	      cufftResult e = cufftExecR2C(plan_e, (float * )d_signal_transformed_e, (cufftComplex * )d_signal_fft_e);
	      
	      if(e != CUFFT_SUCCESS) {
		LOG(log_level::error, "Could not cufftExecR2C in aa_host_acceleration.cu");
	      }
	      
	      power_gpu(event_e, stream_e, samps, a, d_signal_fft_e, d_signal_power_e);
	    }

	  // Threshold odd f-fdot plane
	  cudaStreamSynchronize(stream_o);
	  stats_gpu(event_o, stream_o, samps, &mean, &stddev, h_signal_power_o, d_signal_power_o);
	}

      //Destroy CUFFT context
      cufftDestroy(plan_e);
      cufftDestroy(plan_o);

      //Destroy streams
      gpuErrchk(cudaStreamDestroy(stream_e));
      gpuErrchk(cudaStreamDestroy(stream_o));

      // cleanup even memory
      cudaFreeHost(h_signal_e);
      cudaFreeHost(h_signal_power_e);
      cudaFree(d_signal_in_e);
      cudaFree(d_signal_fft_e);
      cudaFree(d_signal_power_e);
      cudaFree(d_signal_transformed_e);

      // cleanup odd memory
      cudaFreeHost(h_signal_o);
      cudaFreeHost(h_signal_power_o);
      cudaFree(d_signal_in_o);
      cudaFree(d_signal_fft_o);
      cudaFree(d_signal_power_o);
      cudaFree(d_signal_transformed_o);
    }
  }
} //namespace astroaccelerate
