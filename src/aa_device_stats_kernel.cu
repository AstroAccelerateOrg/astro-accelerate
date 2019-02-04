#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "aa_params.hpp"

namespace astroaccelerate {

  /** \brief Sets stats. */
  __global__ void stats_kernel(int half_samps, float *d_sum, float *d_sum_square, float *d_signal_power)
  {

    int t = blockIdx.x * blockDim.x * STATSLOOP + threadIdx.x;

    float local = 0.0;
    float sum = 0.0;
    float sum_square = 0.0;

    for (int i = t; i < t + STATSLOOP * blockDim.x; i += blockDim.x)
      {
	local = d_signal_power[i];
	sum += local;
	sum_square += local * local;
      }
    d_sum[blockIdx.x * blockDim.x + threadIdx.x] = sum;
    d_sum_square[blockIdx.x * blockDim.x + threadIdx.x] = sum_square;
  }

  /** Kernel wrapper function for stats_kernel kernel function. */
  void call_kernel_stats_kernel(const dim3 &block_size, const dim3 &grid_size, const int &smem_bytes, const cudaStream_t &stream,
				const int &half_samps, float *const d_sum, float *const d_sum_square, float *const d_signal_power) {
    stats_kernel<<<block_size, grid_size, smem_bytes, stream>>>(half_samps, d_sum, d_sum_square, d_signal_power);
  }

} //namespace astroaccelerate
