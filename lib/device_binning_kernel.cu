#ifndef BIN_KERNEL_H_
#define BIN_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "AstroAccelerate/params.h"

#define BINARRAYSIZE 2 * BINDIVINT * BINDIVINF

__shared__ float f_line_bin[BINARRAYSIZE];

__global__ void bin(unsigned short *d_input, float *d_output, int in_nsamp)
{
	
	int	c = ((blockIdx.y*BINDIVINF) + threadIdx.y);
	int	out_nsamp = (in_nsamp) / 2;
	int	t_out =  ((blockIdx.x*BINDIVINT) + threadIdx.x);
	int	t_in = 2*t_out;
	
	int	shift_one = ((c*out_nsamp) + t_out);
	int	shift_two = ((c*in_nsamp)  + t_in);

	d_output[(shift_one)] = (float)((d_input[(shift_two)] + d_input[shift_two + 1]));

}

#endif
