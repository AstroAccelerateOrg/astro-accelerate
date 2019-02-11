#include "aa_device_SNR_limited_inplace_kernel.hpp"

namespace astroaccelerate {

  //__device__ __constant__ float c_sqrt_taps[PD_MAXTAPS+1];

  __global__ void PD_ZC_GPU(float *d_input, float *d_output, int maxTaps, int nTimesamples, int nLoops) {
    int x_r, y_r, x_w, y_w;
    int Elements_per_block=PD_NTHREADS*PD_NWINDOWS;
	
    //read
    y_r=(blockIdx.y*blockDim.y + threadIdx.y)*nTimesamples;
    x_r=(blockIdx.x+1)*Elements_per_block + threadIdx.x;
	
    //write
    y_w=(blockIdx.y*blockDim.y + threadIdx.y)*(maxTaps-1)*gridDim.x;
    x_w=blockIdx.x*(maxTaps-1) + threadIdx.x;
	
    for(int f=0; f<nLoops; f++){
      if(x_r<nTimesamples && threadIdx.x<(maxTaps-1)){
	d_output[x_w + y_w + f*WARP]=d_input[x_r + y_r + f*WARP];
      }
    }
  }


  __global__ void PD_GPUv1_const(float *d_input, float *d_temp, unsigned char *d_output_taps, int maxTaps, int nTimesamples, float signal_mean, float signal_sd) {
    extern __shared__ float s_input[]; //dynamically allocated memory for now
	
    int f, gpos_y, gpos_x, spos, itemp;
    float res_SNR[PD_NWINDOWS], SNR, FIR_value, ftemp;
    int res_Taps[PD_NWINDOWS];
	
    //----------------------------------------------
    //----> Reading data
    gpos_y = blockIdx.y*nTimesamples;
    gpos_x = blockIdx.x*PD_NTHREADS*PD_NWINDOWS + threadIdx.x;
    spos=threadIdx.x;
    for(f=0; f<PD_NWINDOWS; f++){
      if( gpos_x<nTimesamples ) {
	s_input[spos]=d_input[gpos_y + gpos_x];
      }
      spos   =  spos  + blockDim.x;
      gpos_x = gpos_x + blockDim.x;
    }
	
    //----> Loading shared data
    itemp=PD_NTHREADS*PD_NWINDOWS + maxTaps -  1;
    gpos_y=blockIdx.y*(maxTaps-1)*gridDim.x;
    gpos_x=blockIdx.x*(maxTaps-1) + threadIdx.x;
    while(spos<itemp){ // && gpos_x<((maxTaps-1)*gridDim.x)
      s_input[spos]=d_temp[gpos_y + gpos_x];
      spos   =  spos  + blockDim.x;
      gpos_x = gpos_x + blockDim.x;
    }
	
    //----> SNR for nTaps=1
    spos=threadIdx.x;
    for(f=0; f<PD_NWINDOWS; f++){
      res_SNR[f]=(s_input[spos]-signal_mean)/signal_sd;
      res_Taps[f]=1;
      spos = spos + blockDim.x;
    }
	
    __syncthreads();
	
    //----------------------------------------------
    //----> FIR calculation loop
    for(f=2; f<=maxTaps; f++){
      ftemp = c_sqrt_taps[f]*signal_sd;
      for(int i=0; i<PD_NWINDOWS; i++){
	spos = threadIdx.x + i*blockDim.x;
	FIR_value=0;
	for(int t=0; t<f; t++){
	  FIR_value+=s_input[spos + t];
	}
	SNR=(FIR_value-f*signal_mean)/(ftemp);		
	if(SNR>res_SNR[i]) {
	  res_SNR[i]=SNR;
	  res_Taps[i]=f;
	}
      }
    }

    //----------------------------------------------
    //---- Writing data
    gpos_y = blockIdx.y*nTimesamples;
    gpos_x = blockIdx.x*PD_NTHREADS*PD_NWINDOWS + threadIdx.x;
    for(int i=0; i<PD_NWINDOWS; i++){
      if( gpos_x<(nTimesamples) ) {
	d_input[gpos_y + gpos_x]=res_SNR[i];
	d_output_taps[gpos_y + gpos_x]=res_Taps[i];
      }
      gpos_x = gpos_x + blockDim.x;
    }
  }

  /** \brief Kernel wrapper function to PD_ZC_GPU kernel function. */
  void call_kernel_PD_ZC_GPU(float *const d_input, float *const d_output, const int &maxTaps, const int &nTimesamples, const int &nLoops) {
    //This kernel looks like it is not called elsewhere in astro-accelerate, hence no implementation is provided.
  }

  /** Kernel wrapper function to PD_GPUv1_const kernel function. */
  void call_kernel_PD_GPUv1_const(float *const d_input, float *const d_temp, unsigned char *const d_output_taps,
				  const int &maxTaps, const int &nTimesamples, const float &signal_mean, const float &signal_sd) {
    //This kernel looks like it is not called elsewhere in astro-accelerate, hence no implementation is provided.
  }


} //namespace astroaccelerate
