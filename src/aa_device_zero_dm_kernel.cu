#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include "aa_device_cuda_deprecated_wrappers.cuh"

namespace astroaccelerate {

  __device__ __inline__ void Reduce_SM(float *s_input){
    for (int i = ( blockDim.x >> 1 ); i > 16; i = i >> 1) {
      if (threadIdx.x < i) {
        s_input[threadIdx.x] = s_input[threadIdx.x] + s_input[threadIdx.x + i];
      }
      __syncthreads();
    }
  }
  
  __device__ __inline__ void Reduce_SM(int *s_input){
    for (int i = ( blockDim.x >> 1 ); i > 16; i = i >> 1) {
      if (threadIdx.x < i) {
        s_input[threadIdx.x] = s_input[threadIdx.x] + s_input[threadIdx.x + i];
      }
      __syncthreads();
    }
  }

  __device__ __inline__ void Reduce_WARP(float *sum){
    for (int q = 16; q > 0; q = q >> 1) {
      *sum += aa_shfl_down(AA_ASSUME_MASK,(*sum), q);
    }
  }
  
  __device__ __inline__ void Reduce_WARP(int *sum){
    for (int q = 16; q > 0; q = q >> 1) {
      *sum += aa_shfl_down(AA_ASSUME_MASK,(*sum), q);
    }
  }

  /**
   * \brief zero_dm_kernel.
   * \todo Needs cleaning and optimizing (WA 21/10/16).
   */
  __global__ void zero_dm_kernel(unsigned short *d_input, int nchans, int nsamp, int nbits, float *normalization_factor) {
    int t  = threadIdx.x;
    extern __shared__ float s_input[];
    float sum = 0.0f;
    int n_iterations = (nchans+blockDim.x-1)/blockDim.x;

    for(int c = 0; c < n_iterations; c++){
      size_t pos = c*blockDim.x + t;
      if (pos < nchans) {
        size_t global_pos = blockIdx.x*((size_t) nchans) + pos;
        sum += d_input[global_pos];
      }
    }

    s_input[t] = sum;
    __syncthreads();

    Reduce_SM(s_input);
    sum = s_input[t];

    Reduce_WARP(&sum);

    if (t == 0) {
      s_input[0] = sum;
    }
    __syncthreads();
    sum = s_input[0];

    sum = sum/((float) nchans);

    for(int c = 0; c < n_iterations; c++){
      size_t pos = c*blockDim.x + t;
      if(pos < nchans) {
        size_t global_pos = blockIdx.x*((size_t) nchans) + pos;
        unsigned short value;
        float result = (float)d_input[global_pos] - sum + normalization_factor[pos];
        if( nbits == 4 ){
          if(result<0) value = 0;
          else if(result>15) value = 15;
          else value = (unsigned short) result;
        }
        else if( nbits == 8 ){
          if(result<0) value = 0;
          else if(result>255) value = 255;
          else value = (unsigned short) result;
        }
        else if( nbits == 16 ){
          if(result<0) value = 0;
          else if(result>65535) value = 65535;
          else value = (unsigned short) result;
        }
        d_input[global_pos] = value;
      }
    }

  }
  
  
  template<typename InputType>
  __global__ void normalization_transposed_kernel(InputType *d_input, size_t nTimesamples, size_t nDMs, float normalization_factor, int nbits){
      __shared__ float s_sum_DDTR_dm[1024];
      int nIterations = (nDMs + 32 - 1)/32;
      float sum = 0;
      size_t pos_t;
      
      // Find partial sum for fix time-sample
      pos_t  = 32*blockIdx.x + threadIdx.x;
      for(int i = 0; i < nIterations; i++){
          size_t pos_dm = 32*i + threadIdx.y;
          if( pos_dm < nDMs && pos_t < nTimesamples){
              size_t global_pos = pos_dm*nTimesamples + pos_t;
              sum = sum + (float) d_input[global_pos];
          }
      }
      
      // Save partial results into shared memory
      s_sum_DDTR_dm[threadIdx.y*32 + threadIdx.x] = sum;
      __syncthreads();
      
      // Find final sum
      for (int i = 16; i > 0; i = i >> 1) {
          int pos = threadIdx.y*32 + threadIdx.x;
          if (threadIdx.y < i) {
              s_sum_DDTR_dm[pos] = s_sum_DDTR_dm[pos] + s_sum_DDTR_dm[pos + i*32];
          }
          __syncthreads();
      }
      
      // Sharing sum value among all threads with same x coordinate
      sum = s_sum_DDTR_dm[threadIdx.x]/((float) nDMs) - normalization_factor;
      
      // Normalization
      pos_t  = 32*blockIdx.x + threadIdx.x;
      for(int i = 0; i < nIterations; i++){
          size_t pos_dm = 32*i + threadIdx.y;
          if( pos_dm < nDMs && pos_t < nTimesamples){
              size_t global_pos = pos_dm*nTimesamples + pos_t;
              InputType value;
              float result  = (float) d_input[global_pos] - sum;
              
              if( nbits == 4 ){
                  if(result<0) value = 0;
                  else if(result>15) value = 15;
                  else value = (InputType) result;
              }
              else if( nbits == 8 ){
                  if(result<0) value = 0;
                  else if(result>255) value = 255;
                  else value = (unsigned short) result;
              }
              else if( nbits == 16 ){
                  if(result<0) value = 0;
                  else if(result>65535) value = 65535;
                  else value = (unsigned short) result;
              }
              else if( nbits == 32 ){
                  value = result;
              }
              
              d_input[global_pos] = value;
          }
      }
  }

  /** \brief Kernel wrapper function for zero_dm_kernel kernel function. */
  void call_kernel_zero_dm_kernel(
    const dim3 &block_size, 
    const dim3 &grid_size, 
    const int &sharedMemory_size,
    unsigned short *const d_input, 
    const int &nchans, 
    const int &nsamp, 
	const int &nbits,
    float *normalization_factor
  ) {
    zero_dm_kernel<<<block_size, grid_size, sharedMemory_size >>>(d_input, nchans, nsamp, nbits, normalization_factor);
  }
  
  void call_kernel_zero_dm_normalization_dm(
      dim3 &nBlocks_per_grid, 
      dim3 &nThreads_per_block, 
      int shared_memory_size, 
      unsigned short *d_input, 
      size_t nTimesamples, 
      size_t nDMs,
      float normalization_factor,
	  int nbits
  ){
    normalization_transposed_kernel<<< nBlocks_per_grid, nThreads_per_block, shared_memory_size>>>(d_input, nTimesamples, nDMs, normalization_factor, nbits);
  }
  
  void call_kernel_post_DDTR_normalization_dm(
      dim3 &nBlocks_per_grid, 
      dim3 &nThreads_per_block, 
      int shared_memory_size, 
      float *d_input, 
      size_t nTimesamples, 
      size_t nDMs,
      float normalization_factor,
	  int nbits
  ){
    normalization_transposed_kernel<<< nBlocks_per_grid, nThreads_per_block, shared_memory_size>>>(d_input, nTimesamples, nDMs, normalization_factor, nbits);
  }
} //namespace astroaccelerate
