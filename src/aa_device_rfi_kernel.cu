#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"


namespace astroaccelerate {

  //{{{ zero dm kernel - needs cleaning and optimizing // WA 21/10/16
  __global__ void rfi_gpu_kernel(unsigned short *d_input, int nchans, int nsamp)
  {

    int c  = blockIdx.x * blockDim.x + threadIdx.x;

    int count =0;

    float stdev = 1000000.0f;
    float mean = 0.0f;
    float sum = 0.0f;
    float sum_squares = 0.0f;
    float cutoff = (4.0f * stdev); 

    for(int out=0; out<4; out++) {
      sum = 0.0f;
      sum_squares = 0.0f;
      count = 0;

      for(int t = 0; t < nsamp; t++) {
	float data=(float)d_input[c*nsamp + t];
	if(data < (mean + cutoff) && data > (mean - cutoff) ) {
	  sum+=data;
	  sum_squares+=(data*data);
	  count++;
	}
      }
      mean = (sum/(float)count);
      sum_squares = ((sum_squares / count) - (mean * mean));
      stdev = sqrt(sum_squares);
      cutoff = (4.0f * stdev); 
    }
	
    for(int t = 0; t < nsamp-4; t++) {
      float data=0.0f;
      for(int x = 0; x<4; x++) data+=(float)d_input[c*nsamp + t + x];
      data=data*0.25f;
      //float data=(float)d_input[c*nsamp + t];
      if(data > (mean + cutoff) || data < (mean - cutoff)) {
	d_input[c*nsamp + t]=(unsigned short)mean;
      }
    }

  }


  /** \brief Kernel wrapper function for rfi_gpu_kernel kernel function. */
  void call_kernel_rfi_gpu_kernel(const dim3 &block_size, const dim3 &grid_size,
				  unsigned short *const d_input, const int &nchans, const int &nsamp) {
    rfi_gpu_kernel<<<block_size, grid_size>>>(d_input, nchans, nsamp);
  }


  //}}}

} //namespace astroaccelerate
