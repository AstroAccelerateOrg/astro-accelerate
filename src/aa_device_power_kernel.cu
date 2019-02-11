#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "aa_params.hpp"

namespace astroaccelerate {

  //{{{ Set stretch
  __global__ void power_kernel(int half_samps, int acc, cufftComplex *d_signal_fft, float *d_signal_power) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;

    if (t < half_samps)
      d_signal_power[t + acc * ( half_samps )] = ( d_signal_fft[t + 1].x * d_signal_fft[t + 1].x + d_signal_fft[t + 1].y * d_signal_fft[t + 1].y );
  }

  //}}}


  __global__ void GPU_simple_power_and_interbin_kernel(float2 *d_input_complex, float *d_output_power, float *d_output_interbinning, int nTimesamples, float norm){
    int pos_x = blockIdx.x*blockDim.x + threadIdx.x;
    int pos_y = blockIdx.y*((nTimesamples>>1)+1);
	
    float2 A, B;
    A.x = 0; A.y = 0; B.x = 0; B.y = 0;
	
    if ( (pos_x < (nTimesamples>>1)) && (pos_x > 0) ) {
      A = d_input_complex[pos_y + pos_x];
      B = d_input_complex[pos_y + pos_x + 1];
		
      A.x = A.x/norm;
      A.y = A.y/norm;
      B.x = B.x/norm;
      B.y = B.y/norm;
    }
	
    if ( (pos_x < (nTimesamples>>1)) ){
      d_output_power[blockIdx.y*(nTimesamples>>1) + pos_x] = A.x*A.x + A.y*A.y;
      d_output_interbinning[blockIdx.y*nTimesamples + 2*pos_x] = A.x*A.x + A.y*A.y;
      d_output_interbinning[blockIdx.y*nTimesamples + 2*pos_x + 1] = 0.616850275f*( (A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y) );
    }
  }

  /** Kernel wrapper function for power_kernel kernel function. */
  void call_kernel_power_kernel(const dim3 &block_size, const dim3 &grid_size, const int &smem_bytes, const cudaStream_t &stream,
				const int &half_samps, const int &acc, cufftComplex *const d_signal_fft, float *const d_signal_power) {
    power_kernel<<<block_size, grid_size, smem_bytes, stream>>>(half_samps, acc, d_signal_fft, d_signal_power);
  }

  /** Kernel wrapper function for GPU_simple_power_and_interbin_kernel kernel function. */
  void call_kernel_GPU_simple_power_and_interbin_kernel(const dim3 &grid_size, const dim3 &block_size,
							float2 *const d_input_complex, float *const d_output_power, float *const d_output_interbinning, const int &nTimesamples, const float &norm) {
    GPU_simple_power_and_interbin_kernel<<<grid_size, block_size>>>(d_input_complex, d_output_power, d_output_interbinning, nTimesamples, norm);
  }

} //namespace astroaccelerate
