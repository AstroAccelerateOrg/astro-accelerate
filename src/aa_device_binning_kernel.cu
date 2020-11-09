#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include <stdio.h>

namespace astroaccelerate {

#define BINARRAYSIZE 2 * BINDIVINT * BINDIVINF

  __shared__ float f_line_bin[BINARRAYSIZE];

  //{{{ bin
  /*
    __global__ void bin(float *bin_buffer, float *input_buffer, int in_nsamp) {

    int	idx = (threadIdx.x + (threadIdx.y * BINDIVINT));
    int	c = ((blockIdx.y * BINDIVINF) + threadIdx.y);

    f_line_bin[idx] = input_buffer[(c*in_nsamp)+(blockIdx.x * BINDIVINT*2)+idx];
    f_line_bin[idx+(BINDIVINT*BINDIVINF)] = input_buffer[(c*in_nsamp)+(blockIdx.x * BINDIVINT*2)+idx+(BINDIVINT*BINDIVINF)];
    __syncthreads();

    int	out_nsamp = in_nsamp / 2;
    int	t_out =  ( (blockIdx.x * BINDIVINT) + threadIdx.x);

    int	shift_one = ((c*out_nsamp) + t_out);
    int	shift_three = (2*threadIdx.x);

    bin_buffer[(shift_one)] = (f_line_bin[(shift_three)] + f_line_bin[shift_three + 1])/2;
    }
  */

  __global__ void bin(unsigned short *d_input, float *d_output, int in_nsamp) {

    int c = ( ( blockIdx.y * BINDIVINF ) + threadIdx.y );
    int out_nsamp = ( in_nsamp ) / 2;
    int t_out = ( ( blockIdx.x * BINDIVINT ) + threadIdx.x );
    int t_in = 2 * t_out;

    size_t shift_one = ( (size_t)(c)*(size_t)(out_nsamp) + (size_t)t_out );
    size_t shift_two = ( (size_t)(c)*(size_t)(in_nsamp)  + (size_t)t_in );

    d_output[( shift_one )] = (float) ( ( d_input[( shift_two )] + d_input[(size_t)(shift_two + 1)] )/2.0f );

  }


  __global__ void DiT_GPU_v2(float const* __restrict__ d_input, float *d_output, unsigned int nDMs, unsigned int nTimesamples, unsigned int dts) {
    float2 ftemp2;
    unsigned int posx, posy, itemp;
	
    posy = (blockIdx.y*DIT_YSTEP + threadIdx.y);
    posx = (blockIdx.x*DIT_ELEMENTS_PER_THREAD*blockDim.x);
	
    //#pragma unroll
    for(int f=0; f<DIT_ELEMENTS_PER_THREAD; f++){
      itemp = (posx + threadIdx.x + f*blockDim.x);
      if( (2*itemp+1)<nTimesamples ){
	ftemp2.x = d_input[posy*nTimesamples + 2*itemp];
	ftemp2.y = d_input[posy*nTimesamples + 2*itemp+1];
	d_output[posy*dts + itemp] = ftemp2.x + ftemp2.y;
      }
    }
  }

  /** \brief Kernel wrapper function for bin kernel function. */
  void call_kernel_bin(const dim3 &num_blocks, const dim3 &threads_per_block, unsigned short *const d_input, float *const d_output, const int &in_nsamp) {
	printf("timesamples: %d\n", in_nsamp);
    bin<<<num_blocks, threads_per_block>>>(d_input, d_output, in_nsamp);
  }

  /** \brief Kernel wrapper function for DiT_GPU_v2 kernel function. */
  void call_kernel_DiT_GPU_v2(const dim3 &gridSize, const dim3 &blockSize, float const *const d_input, float *const d_output, const unsigned int &nDMs, const unsigned int &nTimesamples, const unsigned int &dts) {
    DiT_GPU_v2<<<gridSize,blockSize>>>(d_input, d_output, nDMs, nTimesamples, (nTimesamples>>1));
  }

} //namespace astroaccelerate
  
//}}}
