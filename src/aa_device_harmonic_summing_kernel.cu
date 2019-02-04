#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include "aa_device_harmonic_summing_kernel.hpp"

namespace astroaccelerate {

  __global__ void PHS_GPU_kernel_old(float const* __restrict__ d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics){
    float signal_mean=d_MSD[0];
    float signal_sd=d_MSD[1];
    float HS_value, temp_SNR, SNR;
    ushort max_SNR_harmonic;
    int pos;

    // reading 0th harmonic, i.e. fundamental frequency
    pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
    if( (blockIdx.y*blockDim.x + threadIdx.x)<nSpectra ){
      HS_value = __ldg(&d_input[pos]);
      SNR = (HS_value - signal_mean)/(signal_sd);
      max_SNR_harmonic = 0;
		
      if(blockIdx.x>0) {
	for(int f=1; f<nHarmonics; f++){
	  if( (blockIdx.x + f*blockIdx.x)<nTimesamples ){
	    pos = (blockIdx.x + f*blockIdx.x)*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
	    HS_value = HS_value + __ldg(&d_input[pos]);
	    temp_SNR = __frsqrt_rn(f+1)*(HS_value - signal_mean*(f+1))/(signal_sd); //assuming white noise 
	    if(temp_SNR > SNR){
	      SNR = temp_SNR;
	      max_SNR_harmonic = f;
	    }
	  }
	}
      }
		
      pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
      d_output_SNR[pos] = SNR;
      d_output_harmonics[pos] = max_SNR_harmonic;
    }
  }


  __global__ void PHS_GPU_kernel(float const* __restrict__ d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics){
    //extern __shared__ float s_MSD[];
    float HS_value, temp_SNR, SNR;
    ushort max_SNR_harmonic;
    int pos;

    //int itemp = (2*nHarmonics)/blockDim.x;
    //for(int f=0; f<itemp; f++){
    //	pos = f*blockDim.x + threadIdx.x;
    //	if(pos<nHarmonics) s_MSD[pos] = d_MSD[pos];
    //}

    // reading 0th harmonic, i.e. fundamental frequency
    pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
    if( (blockIdx.y*blockDim.x + threadIdx.x)<nSpectra ){
      HS_value = __ldg(&d_input[pos]);
      SNR = (HS_value - __ldg(&d_MSD[0]))/(__ldg(&d_MSD[1]));
      max_SNR_harmonic = 0;
		
      if(blockIdx.x>0) {
	for(int f=1; f<nHarmonics; f++) {
	  if( (blockIdx.x + f*blockIdx.x)<nTimesamples ) {
	    pos = (blockIdx.x + f*blockIdx.x)*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
	    HS_value = HS_value + __ldg(&d_input[pos]);
	    temp_SNR = (HS_value - __ldg(&d_MSD[f*2]))/(__ldg(&d_MSD[2*f+1])); //assuming white noise 
	    if(temp_SNR > SNR) {
	      SNR = temp_SNR;
	      max_SNR_harmonic = f;
	    }
	  }
	}
      }
		
      pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
      d_output_SNR[pos] = SNR;
      d_output_harmonics[pos] = max_SNR_harmonic;
    }
  }

  /** \brief Kernel wrapper function for PHS_GPU_kernel_old kernel function. */
  void call_kernel_PHS_GPU_kernel_old(const dim3 &grid_size, const dim3 &block_size,
				      float const *const d_input, float *const d_output_SNR, ushort *const d_output_harmonics,
				      float *const d_MSD, const int &nTimesamples, const int &nSpectra, const int &nHarmonics) {
    PHS_GPU_kernel_old<<<grid_size, block_size>>>(d_input, d_output_SNR, d_output_harmonics, d_MSD, nTimesamples, nSpectra, nHarmonics);
  }

  /** \brief Kernel wrapper function for PHS_GPU_kernel kernel function. */
  void call_kernel_PHS_GPU_kernel(const dim3 &grid_size, const dim3 &block_size,
				  float const *const d_input, float *const d_output_SNR, ushort *const d_output_harmonics,
				  float *const d_MSD, const int &nTimesamples, const int &nSpectra, const int &nHarmonics) {
    PHS_GPU_kernel<<<grid_size, block_size>>>(d_input, d_output_SNR, d_output_harmonics, d_MSD, nTimesamples, nSpectra, nHarmonics);
  }

} //namespace astroaccelerate
