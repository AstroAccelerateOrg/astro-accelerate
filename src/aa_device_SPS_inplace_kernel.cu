// Added by Karel Adamek
#include "aa_device_SPS_inplace_kernel.hpp"

namespace astroaccelerate {

  __global__ void PD_ZC_GPU_KERNEL(float *d_input, float *d_output, int maxTaps, int nTimesamples, int nLoops)
  {
    int x_r, y_r, x_w, y_w;
    int Elements_per_block = PD_NTHREADS * PD_NWINDOWS;

    //read
    y_r = ( blockIdx.y * blockDim.y + threadIdx.y ) * nTimesamples;
    x_r = ( blockIdx.x + 1 ) * Elements_per_block + threadIdx.x;

    //write
    y_w = ( blockIdx.y * blockDim.y + threadIdx.y ) * ( maxTaps - 1 ) * gridDim.x;
    x_w = blockIdx.x * ( maxTaps - 1 ) + threadIdx.x;

    for (int f = 0; f < nLoops; f++)
      {
	if (x_r < nTimesamples && threadIdx.x < ( maxTaps - 1 ))
	  {
	    d_output[x_w + y_w + f * WARP] = d_input[x_r + y_r + f * WARP];
	  }
      }
  }

  __global__ void PD_INPLACE_GPU_KERNEL(float *d_input, float *d_temp, unsigned char *d_output_taps, float *d_MSD, int maxTaps, int nTimesamples)
  {
    extern __shared__ float s_input[]; //dynamically allocated memory for now

    int f, i, gpos_y, gpos_x, spos, itemp;
    float res_SNR[PD_NWINDOWS], SNR, temp_FIR_value, FIR_value, ftemp;
    int res_Taps[PD_NWINDOWS];
    float signal_mean, signal_sd, modifier;
    signal_mean = d_MSD[0];
    signal_sd = d_MSD[2];
    modifier = d_MSD[1];

    //----------------------------------------------
    //----> Reading data
    gpos_y = blockIdx.y * nTimesamples;
    gpos_x = blockIdx.x * PD_NTHREADS * PD_NWINDOWS + threadIdx.x;
    spos = threadIdx.x;
    for (f = 0; f < PD_NWINDOWS; f++)
      {
	if (gpos_x < nTimesamples)
	  {
	    s_input[spos] = d_input[gpos_y + gpos_x];
	  }
	spos = spos + blockDim.x;
	gpos_x = gpos_x + blockDim.x;
      }

    //----> Loading shared data
    itemp = PD_NTHREADS * PD_NWINDOWS + maxTaps - 1;
    gpos_y = blockIdx.y * ( maxTaps - 1 ) * gridDim.x;
    gpos_x = blockIdx.x * ( maxTaps - 1 ) + threadIdx.x;
    while (spos < itemp)
      { // && gpos_x<((maxTaps-1)*gridDim.x)
	s_input[spos] = d_temp[gpos_y + gpos_x];
	spos = spos + blockDim.x;
	gpos_x = gpos_x + blockDim.x;
      }

    __syncthreads();

    //----> SNR for nTaps=1
    spos = PD_NWINDOWS * threadIdx.x;
    for (i = 0; i < PD_NWINDOWS; i++)
      {
	res_SNR[i] = ( s_input[spos + i] - signal_mean ) / signal_sd;
	res_Taps[i] = 1;
      }

    //----------------------------------------------
    //----> FIR calculation loop
    FIR_value = s_input[spos];
    for (f = 1; f < maxTaps; f++)
      {
	//nTaps=f+1;!
	ftemp = signal_sd + f * modifier;
	spos = PD_NWINDOWS * threadIdx.x;

	// 0th element from NWINDOW
	i = 0;
	FIR_value += s_input[spos + f];

	SNR = ( FIR_value - ( f + 1 ) * signal_mean ) / ( ftemp );
	if (SNR > res_SNR[i])
	  {
	    res_SNR[i] = SNR;
	    res_Taps[i] = f + 1;
	  }

	temp_FIR_value = FIR_value;
	for (i = 1; i < PD_NWINDOWS; i++)
	  {
	    temp_FIR_value = temp_FIR_value - s_input[spos + i - 1] + s_input[spos + f + i];

	    SNR = ( temp_FIR_value - ( f + 1 ) * signal_mean ) / ( ftemp );
	    if (SNR > res_SNR[i])
	      {
		res_SNR[i] = SNR;
		res_Taps[i] = f + 1;
	      }
	  }
      }

    //----------------------------------------------
    //---- Writing data
    gpos_y = blockIdx.y * nTimesamples;
    gpos_x = blockIdx.x * PD_NTHREADS * PD_NWINDOWS + PD_NWINDOWS * threadIdx.x;
    for (i = 0; i < PD_NWINDOWS; i++)
      {
	if (( gpos_x + i ) < ( nTimesamples ))
	  {
	    d_input[gpos_y + gpos_x + i] = res_SNR[i];
	    d_output_taps[gpos_y + gpos_x + i] = res_Taps[i];
	  }
      }
  }

  /** \brief Kernel wrapper function to PD_ZC_GPU_KERNEL kernel function. */
  void call_kernel_PD_ZC_GPU_KERNEL(const dim3 &grid_size, const dim3 &block_size, float *const d_input, float *const d_output, const int &maxTaps, const int &nTimesamples, const int &nLoops) {
    PD_ZC_GPU_KERNEL<<<grid_size, block_size>>>(d_input, d_output, maxTaps, nTimesamples, nLoops);
  }

  /** \brief Kernel wrapper function to PD_INPLACE_GPU_KERNEL kernel function. */
  void call_kernel_PD_INPLACE_GPU_KERNEL(const dim3 &grid_size, const dim3 &block_size, const int &SM_size, float *const d_input,
					 float *const d_temp, unsigned char *const d_output_taps, float *const d_MSD,
					 const int &maxTaps, const int &nTimesamples) {
    PD_INPLACE_GPU_KERNEL<<<grid_size, block_size, SM_size>>>(d_input, d_temp, d_output_taps,
							      d_MSD, maxTaps, nTimesamples);
  }

} //namespace astroaccelerate
