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

  __device__ __inline__ void Reduce_WARP(float *sum){
    for (int q = 16; q > 0; q = q >> 1) {
      *sum += aa_shfl_down(AA_ASSUME_MASK,(*sum), q);
    }
  }

  /**
   * \brief zero_dm_kernel.
   * \todo Needs cleaning and optimizing (WA 21/10/16).
   */
  __global__ void zero_dm_kernel(unsigned short *d_input, int nchans, int nsamp, int nbits, float normalization_factor)
  {

    int t  = threadIdx.x;

    extern __shared__ float s_input[];
    float sum = 0.0f;

    int n_iterations = (nchans+blockDim.x-1)/blockDim.x;

    for(int c = 0; c < n_iterations; c++){
      if ((c*blockDim.x + t) < nchans) {
        sum += d_input[blockIdx.x*nchans + c*blockDim.x + t];
      }
    }

    s_input[t] = sum;
    __syncthreads();

    Reduce_SM(s_input);
    sum = s_input[t];

    Reduce_WARP(&sum);

    if (t  == 0) s_input[0] = sum;
    __syncthreads();
    sum = s_input[0];

    sum = (sum/(float)nchans-normalization_factor);

    for(int c = 0; c < n_iterations; c++){
      if ((c*blockDim.x + t) < nchans) {
		unsigned short value;
		float result = (float)d_input[blockIdx.x*nchans + c*blockDim.x + t] - sum;
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
        d_input[blockIdx.x*nchans + c*blockDim.x + t] = value;
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
    const float &normalization_factor
  ) {
    zero_dm_kernel<<<block_size, grid_size, sharedMemory_size >>>(d_input, nchans, nsamp, nbits, normalization_factor);
  }
} //namespace astroaccelerate
