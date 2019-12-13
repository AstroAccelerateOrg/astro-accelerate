#include "aa_params.hpp"
#include "aa_device_MSD_Configuration.hpp"
#include "aa_device_MSD_shared_kernel_functions.cuh"
#include "aa_device_MSD_normal_kernel.hpp"

namespace astroaccelerate {

  // Computes partials for mean and standard deviation of the data with offset at the end
  // PD_THREADS could be replaced it is not required to be #defined
  __global__ void MSD_GPU_limited(float const* __restrict__ d_input, float *d_output, int y_steps, int nTimesamples, int offset) {
    __shared__ float s_input[3*PD_NTHREADS];
    float M, S, j, ftemp;
  
    int spos = blockIdx.x*PD_NTHREADS + threadIdx.x;
    size_t gpos = (size_t)(blockIdx.y*y_steps*nTimesamples) + (size_t)spos;
    M=0;	S=0;	j=0;
    if( spos < (nTimesamples - offset) ){
    
      ftemp=__ldg(&d_input[gpos]);
      Initiate( &M, &S, &j, ftemp);
    
      gpos = gpos + (size_t)nTimesamples;
      for (int yf = 1; yf < y_steps; yf++) {
	ftemp=__ldg(&d_input[gpos]);
	Add_one( &M, &S, &j, ftemp);
	gpos = gpos + (size_t)nTimesamples;
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
      d_output[MSD_PARTIAL_SIZE*gpos] = M;
      d_output[MSD_PARTIAL_SIZE*gpos + 1] = S;
      d_output[MSD_PARTIAL_SIZE*gpos + 2] = j;
    }
  }

  /** \brief Wrapper function to kernel function MSD_GPU_limited. */
  void call_kernel_MSD_GPU_limited(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, const int &y_steps, const int &nTimesamples, const int &offset) {
    MSD_GPU_limited<<<grid_size, block_size>>>(d_input, d_output, y_steps, nTimesamples, offset);
  }

} //namespace astroaccelerate
