#include <cuda.h>
#include <cuda_runtime.h>
#include "params.hpp"
#include "device_set_stretch_kernel.hpp"

//{{{ Set stretch
__global__ void set_stretch_kernel(int samps, float mean, float *d_input) {

	int t = blockIdx.x * blockDim.x + threadIdx.x;

	if (t >= 0 && t < samps)
		d_input[t] = mean;
}

void call_kernel_set_stretch_kernel(dim3 block_size, dim3 grid_size,
				    int smem_bytes, cudaStream_t stream,
				    int samps, float mean, float *d_input) {
  set_stretch_kernel<<<block_size, grid_size, smem_bytes, stream>>>(samps, mean, d_input);
}

//}}}
