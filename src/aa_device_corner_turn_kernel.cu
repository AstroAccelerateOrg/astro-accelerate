#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include <stdio.h>

namespace astroaccelerate {

  __global__ void simple_corner_turn_kernel(float *d_input, float *d_output, int primary_size, int secondary_size){

    size_t primary = blockIdx.x * blockDim.x + threadIdx.x;
    size_t secondary = blockIdx.y * blockDim.y + threadIdx.y;

    d_output[(size_t)primary*secondary_size + secondary] = (float) __ldg(&d_input[(size_t)secondary*primary_size + primary]);
  }


  __global__ void corner_turn_SM_kernel(float const* __restrict__ d_input, float *d_output, int primary_size, int secondary_size) {
    __shared__ float s_input[WARP*(WARP+1)*CT_CORNER_BLOCKS];
	
    int i, spos, itemp, pc, sc;
    size_t gpos;
	
    int warp_id = threadIdx.x>>5;
    int local_id = threadIdx.x & (WARP - 1);
	
    gpos=(size_t)((size_t)(blockIdx.y*(blockDim.x>>5)) + (size_t)warp_id)*CT_ROWS_PER_WARP*primary_size + (size_t)(blockIdx.x*CT_CORNER_BLOCKS*WARP) + (size_t)local_id;
    for(int by=0; by<CT_ROWS_PER_WARP; by++){
      spos=local_id*WARP + local_id + warp_id*CT_ROWS_PER_WARP + by;
      for(int bx=0; bx<CT_CORNER_BLOCKS; bx++){ // temporary 
	s_input[spos]=d_input[gpos];
	gpos=gpos + (size_t)WARP;
	spos=spos + WARP*(WARP+1);
      }
      gpos=gpos + (size_t)primary_size - (size_t)(CT_CORNER_BLOCKS*WARP);
    }
	
    __syncthreads();
	
    itemp=warp_id*CT_ROWS_PER_WARP*CT_CORNER_BLOCKS;
    for(i=0; i<CT_ROWS_PER_WARP*CT_CORNER_BLOCKS; i++){
      pc = (blockIdx.x*CT_CORNER_BLOCKS*WARP + itemp + i);
      sc = WARP*blockIdx.y + local_id;
      if( pc<primary_size && sc<secondary_size ) {
	gpos=(size_t)(pc*secondary_size) + (size_t)sc;
	spos=(itemp + i)*(WARP+1) + local_id;
	d_output[gpos]=s_input[spos];
      }
    }
  }

  __global__ void simple_corner_turn_kernel(unsigned short *d_input, float *d_output, int nchans, int nsamp) {

    size_t t = blockIdx.x * blockDim.x + threadIdx.x;
    size_t c = blockIdx.y * blockDim.y + threadIdx.y;

    d_output[(size_t)(c * nsamp) + t] = (float) __ldg(&d_input[(size_t)(t * nchans) + c]);

  }

  __global__ void swap(unsigned short *d_input, float *d_output, int nchans, int nsamp) {

    size_t t = blockIdx.x * blockDim.x + threadIdx.x;
    size_t c = blockIdx.y * blockDim.y + threadIdx.y;

    d_input[(size_t)(c * nsamp) + t] = (unsigned short) __ldg(&d_output[(size_t)(c * nsamp) + t]);

  }

  /** \brief Kernel wrapper function for simple_corner_turn_kernel kernel function. */
  void call_kernel_simple_corner_turn_kernel(const dim3 &block_size, const dim3 &grid_size, unsigned short *const d_input, float *const d_output, const int &primary_size, const int &secondary_size) {
    simple_corner_turn_kernel<<<block_size, grid_size>>>(d_input, d_output, primary_size, secondary_size);
  }

  /** \brief Kernel wrapper function for simple_corner_turn_kernel kernel function. */
  void call_kernel_simple_corner_turn_kernel(const dim3 &block_size, const dim3 &grid_size, float *const d_input, float *const d_output, const int &primary_size, const int &secondary_size) {
    simple_corner_turn_kernel<<<block_size, grid_size>>>(d_input, d_output, primary_size, secondary_size);
  }

  /** \brief Kernel wrapper function for corner_turn_SM_kernel kernel function. */
  void call_kernel_corner_turn_SM_kernel(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, const int &primary_size, const int &secondary_size) {
    corner_turn_SM_kernel<<<grid_size,block_size>>>(d_input, d_output, primary_size, secondary_size);
  }

  /** \brief Kernel wrapper function for swap kernel function. */
  void call_kernel_swap(const dim3 &block_size, const dim3 &grid_size, unsigned short *const d_input, float *const d_output, const int &nchans, const int &nsamp) {
    swap<<<block_size, grid_size>>>(d_input, d_output, nchans, nsamp);
  }

} //namespace astroaccelerate
