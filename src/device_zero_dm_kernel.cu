#include "params.hpp"
#include <cuda.h>
#include <cuda_runtime.h>

//{{{ zero dm kernel - needs cleaning and optimizing // WA 21/10/16
__global__ void zero_dm_kernel(unsigned short* d_input,
                               int             nchans,
                               int             nsamp,
                               float           normalization_factor) {

  int t = blockIdx.x * blockDim.x + threadIdx.x;

  float sum = 0.0f;
  for(int c = 0; c < nchans; c++)
    sum += (float)__ldg(&d_input[t * nchans + c]);
  sum = (sum / (float)nchans - normalization_factor);
  for(int c = 0; c < nchans; c++)
    d_input[t * nchans + c] =
        (unsigned short)((unsigned char)((float)d_input[t * nchans + c] - sum));
}

void call_kernel_zero_dm_kernel(const dim3&           block_size,
                                const dim3&           grid_size,
                                unsigned short* const d_input,
                                const int&            nchans,
                                const int&            nsamp,
                                const float&          normalization_factor) {
  zero_dm_kernel<<<block_size, grid_size>>>(
      d_input, nchans, nsamp, normalization_factor);
}

//}}}
