#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include <stdio.h>

namespace astroaccelerate {


  template<typename inType>
  __inline__ __device__ void LoadMatrix(inType *s_submatrix, inType const* __restrict__ d_input, size_t primary_size, size_t secondary_size, size_t block_pos_x, size_t block_pos_y){
      size_t warp_id = threadIdx.x>>5;
      size_t local_id = threadIdx.x & (WARP - 1);
      size_t gpos, total_size, spos;
      
      gpos = (block_pos_y*(blockDim.x>>5) + warp_id)*CT_ROWS_PER_WARP*primary_size + block_pos_x*WARP + local_id;
      total_size = primary_size*secondary_size;
      for(int by=0; by<CT_ROWS_PER_WARP; by++){
          spos=local_id*WARP + local_id + warp_id*CT_ROWS_PER_WARP + by;
          
          if(gpos<total_size) s_submatrix[spos] = d_input[gpos];
          else s_submatrix[spos] = 0;
          
          gpos = gpos + primary_size;
      }
  }
  
  template<typename inType>
  __inline__ __device__ void WriteMatrix(inType *d_output, inType *s_submatrix, size_t primary_size, size_t secondary_size, size_t block_pos_x, size_t block_pos_y){
      size_t local_id = threadIdx.x & (WARP - 1);
      size_t itemp = (threadIdx.x>>5)*CT_ROWS_PER_WARP;
      
      for(int i=0; i<CT_ROWS_PER_WARP; i++){
          size_t pc = (block_pos_x*WARP + itemp + i);
          size_t sc = WARP*block_pos_y + local_id;
          if( pc<primary_size && sc<secondary_size ) {
              size_t gpos = pc*secondary_size + sc;
              int spos = (itemp + i)*(WARP+1) + local_id;
              d_output[gpos] = s_submatrix[spos];
          }
      }
  }
  
  template<typename inType>
  __global__ void corner_turn_SM_kernel(inType const* __restrict__ d_input, inType *d_output, size_t primary_size, size_t secondary_size) {
      __shared__ inType s_input[WARP*(WARP+1)*CT_CORNER_BLOCKS];
      
      for(size_t bx=0; bx<CT_CORNER_BLOCKS; bx++){
          LoadMatrix(&s_input[WARP*(WARP+1)*bx], d_input, primary_size, secondary_size, CT_CORNER_BLOCKS*blockIdx.x + bx, blockIdx.y);
      }
      
      __syncthreads();
      
      for(size_t bx=0; bx<CT_CORNER_BLOCKS; bx++){
          WriteMatrix(d_output, &s_input[WARP*(WARP+1)*bx], primary_size, secondary_size, CT_CORNER_BLOCKS*blockIdx.x + bx, blockIdx.y);
      }
  }
  
  
  template<typename inType>
  __global__ void swap_content(inType *d_destination, inType *d_source, size_t primary_size, size_t secondary_size){
      size_t pos_x = blockIdx.x*blockDim.x + threadIdx.x;
      
      for(int f=0; f<WARP; f++) {
          size_t pos_y = blockIdx.y*WARP + f;
          if(pos_y<secondary_size && pos_x<primary_size){
              size_t gpos = pos_y*primary_size + pos_x;
              d_destination[gpos] = d_source[gpos];
          }
      }
  }



  __global__ void simple_corner_turn_kernel(float *d_input, float *d_output, int primary_size, int secondary_size){

    size_t primary = blockIdx.x * blockDim.x + threadIdx.x;
    size_t secondary = blockIdx.y * blockDim.y + threadIdx.y;

    d_output[(size_t)primary*secondary_size + secondary] = (float) __ldg(&d_input[(size_t)secondary*primary_size + primary]);
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

  //-----------------------------------------------------------------------------

  /** \brief Kernel wrapper function for simple_corner_turn_kernel kernel function. */
  void call_kernel_simple_corner_turn_kernel(const dim3 &block_size, const dim3 &grid_size, unsigned short *const d_input, float *const d_output, const int &primary_size, const int &secondary_size) {
    simple_corner_turn_kernel<<<block_size, grid_size>>>(d_input, d_output, primary_size, secondary_size);
  }

  /** \brief Kernel wrapper function for simple_corner_turn_kernel kernel function. */
  void call_kernel_simple_corner_turn_kernel(const dim3 &block_size, const dim3 &grid_size, float *const d_input, float *const d_output, const int &primary_size, const int &secondary_size) {
    simple_corner_turn_kernel<<<block_size, grid_size>>>(d_input, d_output, primary_size, secondary_size);
  }

  /** \brief Kernel wrapper function for corner_turn_SM_kernel kernel function. */
  void call_kernel_corner_turn_SM_kernel(
      const dim3 &grid_size, 
      const dim3 &block_size, 
      float const *const d_input, 
      float *const d_output, 
      const size_t &primary_size, 
      const size_t &secondary_size
  ) {
    corner_turn_SM_kernel<<<grid_size,block_size>>>(d_input, d_output, primary_size, secondary_size);
  }
  
  /** \brief Kernel wrapper function for corner_turn_SM_kernel kernel function. */
  void call_kernel_corner_turn_SM_kernel(
      const dim3 &grid_size, 
      const dim3 &block_size, 
      unsigned short const *const d_input, 
      unsigned short *const d_output, 
      const size_t &primary_size, 
      const size_t &secondary_size
  ) {
    corner_turn_SM_kernel<<<grid_size,block_size>>>(d_input, d_output, primary_size, secondary_size);
  }
  
  void call_kernel_swap_content(
      const dim3 &grid_size, 
      const dim3 &block_size, 
      float *d_destination, 
      float *d_source, 
      size_t primary_size, 
      size_t secondary_size
  ) {
      swap_content<<<grid_size,block_size>>>(d_destination, d_source, primary_size, secondary_size);
  }
  
  void call_kernel_swap_content(
      const dim3 &grid_size, 
      const dim3 &block_size, 
      unsigned short *d_destination, 
      unsigned short *d_source, 
      size_t primary_size, 
      size_t secondary_size
  ) {
      swap_content<<<grid_size,block_size>>>(d_destination, d_source, primary_size, secondary_size);
  }

  /** \brief Kernel wrapper function for swap kernel function. */
  void call_kernel_swap(const dim3 &block_size, const dim3 &grid_size, unsigned short *const d_input, float *const d_output, const int &nchans, const int &nsamp) {
    swap<<<block_size, grid_size>>>(d_input, d_output, nchans, nsamp);
  }

} //namespace astroaccelerate
