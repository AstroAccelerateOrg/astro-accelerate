#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include "aa_device_single_pulse_search_kernel.hpp"
#include "aa_device_SPS_inplace_kernel.hpp"

namespace astroaccelerate {

  __global__ void PD_SEARCH_GPU(float const* __restrict__ d_input, float *d_output, float *d_output_taps, float *d_MSD, int maxTaps, int nTimesamples)
  {
    extern __shared__ float s_input[]; //dynamically allocated memory for now

    int f, spos, itemp;
    size_t gpos;
    float res_SNR[PD_NWINDOWS], SNR, FIR_value, ftemp, signal_mean, signal_sd;
    int res_Taps[PD_NWINDOWS];
    signal_mean = d_MSD[0];
    signal_sd = d_MSD[1];

    //----------------------------------------------
    //----> Reading data
    gpos = (size_t)(blockIdx.y*nTimesamples) + (size_t)(blockIdx.x*PD_NTHREADS*PD_NWINDOWS);
    spos = threadIdx.x;
    for (f = 0; f < PD_NWINDOWS; f++)
      {
	s_input[spos] = __ldg(&d_input[gpos + (size_t)spos]);
	spos = spos + blockDim.x;
      }

    itemp = PD_NTHREADS*PD_NWINDOWS + maxTaps - 1;
    while (spos < itemp)
      {
	s_input[spos] = __ldg(&d_input[gpos + (size_t)spos]);
	spos = spos + blockDim.x;
      }

    //----> SNR for nTaps=1
    spos = threadIdx.x;
    for (f = 0; f < PD_NWINDOWS; f++)
      {
	res_SNR[f] = ( s_input[spos] - signal_mean )/signal_sd;
	res_Taps[f] = 1;
	spos = spos + blockDim.x;
      }

    __syncthreads();

    //----------------------------------------------
    //----> FIR calculation loop
    for (f = 2; f <= maxTaps; f++)
      {
	ftemp = c_sqrt_taps[f]*signal_sd;
	for (int i = 0; i < PD_NWINDOWS; i++)
	  {
	    spos = threadIdx.x + i*blockDim.x;
	    FIR_value = 0;
	    for (int t = 0; t < f; t++)
	      {
		FIR_value += s_input[spos + t];
	      }
	    SNR = ( FIR_value - f*signal_mean )/( ftemp );
	    if (SNR > res_SNR[i])
	      {
		res_SNR[i] = SNR;
		res_Taps[i] = f;
	      }
	  }
      }

    //----------------------------------------------
    //---- Writing data
    gpos = (size_t)(blockIdx.y*nTimesamples) + (size_t)(blockIdx.x*PD_NTHREADS*PD_NWINDOWS);
    spos = threadIdx.x;
    for (int i = 0; i < PD_NWINDOWS; i++)
      {
	d_output[gpos + (size_t)spos] = res_SNR[i];
	d_output_taps[gpos + (size_t)spos] = res_Taps[i];
	spos = spos + blockDim.x;
      }
  }

  /** \brief Kernel wrapper function for PD_SEARCH_GPU kernel function. */
  void call_kernel_PD_SEARCH_GPU(const dim3 &grid_size, const dim3 &block_size, const int &sm_size,
				 float const *const d_input, float *const d_output, float *const d_output_taps, float *const d_MSD, const int &maxTaps, const int &nTimesamples) {
    PD_SEARCH_GPU<<<grid_size, block_size, sm_size>>>(d_input, d_output, d_output_taps, d_MSD, maxTaps, nTimesamples);
  }

} //namespace astroaccelerate
