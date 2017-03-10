#ifndef BIN_KERNEL_H_
#define BIN_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

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

__global__ void bin(unsigned short *d_input, float *d_output, int in_nsamp)
{

	int c = ( ( blockIdx.y * BINDIVINF ) + threadIdx.y );
	int out_nsamp = ( in_nsamp ) / 2;
	int t_out = ( ( blockIdx.x * BINDIVINT ) + threadIdx.x );
	int t_in = 2 * t_out;

	int shift_one = ( ( c * out_nsamp ) + t_out );
	int shift_two = ( ( c * in_nsamp ) + t_in );

	d_output[( shift_one )] = (float) ( ( d_input[( shift_two )] + d_input[shift_two + 1] )/2.0f );

}

#endif

//}}}
