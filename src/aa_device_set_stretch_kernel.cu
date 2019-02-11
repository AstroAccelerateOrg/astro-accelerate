#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include "aa_device_set_stretch_kernel.hpp"

namespace astroaccelerate {
  /** \brief Set stretch kernel. */
  __global__ void set_stretch_kernel(int samps, float mean, float *d_input) {

    int t = blockIdx.x * blockDim.x + threadIdx.x;

    if (t >= 0 && t < samps)
      d_input[t] = mean;
  }

  /** \brief Kernel wrapper function for set_stretch_kernel kernel function. */
  void call_kernel_set_stretch_kernel(const dim3 &block_size, const dim3 &grid_size,
				      const int &smem_bytes, const cudaStream_t &stream,
				      const int &samps, const float &mean, float *const d_input) {
    set_stretch_kernel<<<block_size, grid_size, smem_bytes, stream>>>(samps, mean, d_input);
  }
} //namespace astroaccelerate
