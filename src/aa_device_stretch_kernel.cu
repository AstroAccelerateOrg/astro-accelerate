#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include <cufft.h>

namespace astroaccelerate {
  /** \brief Stretch. */
  __global__ void stretch_kernel(int acc, int samps, float tsamp, float *d_input, float *d_output, float t_zero, float multiplier, float tsamp_inverse) {

    int t = blockIdx.x * blockDim.x + threadIdx.x;

    float p_time = t * ( t_zero + ( multiplier * ( t - 1.0f ) ) );

    int stretch_index = __float2int_rz(p_time * tsamp_inverse);

    if (stretch_index >= 0 && stretch_index < samps)
      d_output[stretch_index] = d_input[t];
  }

  void call_kernel_stretch_kernel(const dim3 &block_size, const dim3 &grid_size, const int &smem_bytes, const cudaStream_t &stream,
				  const int &acc, const int &samps, const float &tsamp, float *const d_input, float *const d_output, const float &t_zero,
				  const float &multiplier, const float &tsamp_inverse) {
    stretch_kernel<<<block_size, grid_size, smem_bytes, stream>>>(acc, samps, tsamp, d_input, d_output, t_zero, multiplier,
								  tsamp_inverse);
  }
} //namespace astroaccelerate
