#ifndef ASTRO_ACCELERATE_DEVICE_MSD_NORMAL_KERNEL_CU
#define ASTRO_ACCELERATE_DEVICE_MSD_NORMAL_KERNEL_CU

#include "params.hpp"
#include "device_MSD_Configuration.hpp"
#include "device_MSD_shared_kernel_functions.cuh"
#include "device_MSD_normal_kernel.hpp"



// Computes partials for mean and standard deviation of the data with offset at the end
// PD_THREADS could be replaced it is not required to be #defined
__global__ void MSD_GPU_limited(float const* __restrict__ d_input, float *d_output, int y_steps, int nTimesamples, int offset) {
  __shared__ float s_input[3*PD_NTHREADS];
  float M, S, j, ftemp;
  
  int spos = blockIdx.x*PD_NTHREADS + threadIdx.x;
  int gpos = blockIdx.y*y_steps*nTimesamples + spos;
  M=0;	S=0;	j=0;
  if( spos<(nTimesamples-offset) ){
    
    ftemp=__ldg(&d_input[gpos]);
    Initiate( &M, &S, &j, ftemp);
    
    gpos = gpos + nTimesamples;
    for (int yf = 1; yf < y_steps; yf++) {
      ftemp=__ldg(&d_input[gpos]);
      Add_one( &M, &S, &j, ftemp);
      gpos = gpos + nTimesamples;
    }
  }
  
  s_input[threadIdx.x] = M;
  s_input[blockDim.x + threadIdx.x] = S;
  s_input[2*blockDim.x + threadIdx.x] = j;
  
  __syncthreads();
  
  Reduce_SM( &M, &S, &j, s_input );
  Reduce_WARP( &M, &S, &j);
  
  //----------------------------------------------
  //---- Writing data
  if (threadIdx.x == 0) {
    gpos = blockIdx.y*gridDim.x + blockIdx.x;
    d_output[3*gpos] = M;
    d_output[3*gpos + 1] = S;
    d_output[3*gpos + 2] = j;
  }
}

void call_kernel_MSD_GPU_limited(dim3 grid_size, dim3 block_size, float const* d_input, float* d_output, int y_steps, int nTimesamples, int offset) {
  MSD_GPU_limited<<<grid_size, block_size>>>(d_input, d_output, y_steps, nTimesamples, offset);
}

#endif
