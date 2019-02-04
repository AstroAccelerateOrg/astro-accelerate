// Added by Karel Adamek 
#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"

namespace astroaccelerate {

  __global__ void PD_FIR_GPU(float const* __restrict__ d_input, float *d_output, int nTaps, int nLoops, int nTimesamples) {
    extern __shared__ float s_input[];

    int itemp, pos;
    float sum[PD_FIR_NWINDOWS];

    //----------------------------------------------
    //---- Reading data
    itemp = PD_FIR_ACTIVE_WARPS * WARP * PD_FIR_NWINDOWS + nTaps - 1;
    for (int i = 0; i < nLoops; i++)
      {
	pos = i * PD_FIR_ACTIVE_WARPS * WARP + threadIdx.x;
	if (pos < itemp)
	  {
	    s_input[pos] = d_input[blockIdx.y * nTimesamples + blockIdx.x * PD_FIR_ACTIVE_WARPS * WARP * PD_FIR_NWINDOWS + pos];
	  }
      }

    __syncthreads();

    //----------------------------------------------
    //---- Calculating FIR version 2

    pos = PD_FIR_NWINDOWS * threadIdx.x;
    sum[0] = 0;
    for (int t = 0; t < nTaps; t++)
      {
	sum[0] += s_input[pos + t];
      }
    for (int i = 1; i < PD_FIR_NWINDOWS; i++)
      {
	pos = PD_FIR_NWINDOWS * threadIdx.x + i - 1;
	sum[i] = sum[i - 1] - s_input[pos] + s_input[pos + nTaps];
      }

    //----------------------------------------------
    //---- Writing data	
    for (int i = 0; i < PD_FIR_NWINDOWS; i++)
      {
	pos = PD_FIR_NWINDOWS * threadIdx.x + i;
	d_output[blockIdx.y * nTimesamples + blockIdx.x * PD_FIR_ACTIVE_WARPS * WARP * PD_FIR_NWINDOWS + pos] = sum[i];
      }
  }



  __global__ void PD_FIR_GPUv1(float const* __restrict__ d_input, float *d_output, int nTaps, int nLoops, unsigned int nTimesamples) {
    extern __shared__ float s_input[];

    unsigned int itemp,pos;
    float sum[PD_FIR_NWINDOWS];


    //----------------------------------------------
    //---- Reading data
    itemp=PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + nTaps - 1; //amount of data needed probably define it in params.h as it is also size of the shared memory
    for(int i=0; i<nLoops; i++){
      pos=i*PD_FIR_ACTIVE_WARPS*WARP + threadIdx.x;
      if(pos<itemp){
	s_input[pos]=d_input[blockIdx.y*nTimesamples + blockIdx.x*PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + pos];
      }
    }
	
    __syncthreads();
	
	
    //----------------------------------------------
    //---- Calculating FIR version 1
    for(int i=0; i<PD_FIR_NWINDOWS; i++){
      pos=i*PD_FIR_ACTIVE_WARPS*WARP + threadIdx.x;
      sum[i]=0;
      for(int t=0; t<nTaps; t++){
	sum[i]+=s_input[pos + t];
      }
    }

    //----------------------------------------------
    //---- Writing data	
    for(int i=0; i<PD_FIR_NWINDOWS; i++){
      pos=i*PD_FIR_ACTIVE_WARPS*WARP + threadIdx.x;
      d_output[blockIdx.y*nTimesamples + blockIdx.x*PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + pos]=sum[i];
    }
    // This so far ignores DM_PER_BLOCK!!!
  }



  __global__ void Fir_L1(float const* __restrict__ d_data, float* d_output, int nTaps, int nTimesamples) {
    int t = 0;
    int posy = blockIdx.y;
    int posx = blockIdx.x*blockDim.x + threadIdx.x;
    float ftemp;

    if(posx<nTimesamples-nTaps+1){
      ftemp=0;
      for(t=0; t<nTaps; t++){
	ftemp += __ldg(&d_data[posy*nTimesamples + posx + t]);
      }

      d_output[posy*nTimesamples + posx] = ftemp;
    }
  }

  /** \brief Kernel wrapper function for PD_FIR_GPU kernel function. */
  void call_kernel_PD_FIR_GPU(const dim3 &grid_size, const dim3 &block_size, const int &SM_size, float const *const d_input, float *const d_output, const int &nTaps, const int &nLoops, const int &nTimesamples) {
    PD_FIR_GPU<<<grid_size, block_size, SM_size>>>(d_input, d_output, nTaps, nLoops, nTimesamples);
  }

  /** \brief Kernel wrapper function for PD_FIR_GPUv1 kernel function. */
  void call_kernel_PD_FIR_GPUv1(const dim3 &grid_size, const dim3 &block_size, const int &SM_size, float const *const d_input, float *const d_output, const int &nTaps, const int &nLoops, const unsigned int &nTimesamples) {
    PD_FIR_GPUv1<<<grid_size, block_size, SM_size>>>(d_input, d_output, nTaps, nLoops, nTimesamples);
  }

  /** \brief Kernel wrapper function for Fir_L1 kernel function. */
  void call_kernel_Fir_L1(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, const int &nTaps, const int &nTimesamples) {
    Fir_L1<<<grid_size, block_size>>>(d_input, d_output, nTaps, nTimesamples);
  }

} //namespace astroaccelerate
