#include "aa_device_SNR_limited_kernel.hpp"

namespace astroaccelerate {

  __global__ void SNR_GPU_limited(float *d_FIR_input, float *d_SNR_output, ushort *d_SNR_taps, float *d_MSD, int x_steps, int nTaps, int nColumns, int offset) {
    int local_id = threadIdx.x & ( WARP - 1 );
    int warp_id = threadIdx.x >> 5;
    int dim_y = blockDim.x >> 5;

    int pos_x, pos_y;
    float old_SNR, new_SNR;

    float signal_mean = d_MSD[0];
    float signal_sd = d_MSD[1];

    pos_y = ( blockIdx.y * dim_y + warp_id ) * nColumns;
    pos_x = blockIdx.x * WARP * x_steps + local_id;

    for (int xf = 0; xf < x_steps; xf++)
      {
	if (pos_x < ( nColumns - offset ))
	  {
	    old_SNR = d_SNR_output[pos_y + pos_x];
	    //new_SNR = (d_FIR_input[pos_y + pos_x]-nTaps*signal_mean)/(sqrt((float) nTaps)*signal_sd);
	    new_SNR = ( d_FIR_input[pos_y + pos_x] - signal_mean ) / ( signal_sd );
	    if (nTaps == 1)
	      {
		//if(new_SNR>3.0){
		//	d_FIR_input[pos_y + pos_x]=d_MSD[0];
		//	d_SNR_output[pos_y + pos_x]=0;
		//	d_SNR_taps[pos_y + pos_x]=nTaps;
		//}
		//else {
		d_SNR_output[pos_y + pos_x] = new_SNR;
		d_SNR_taps[pos_y + pos_x] = nTaps;
		//}
	      }
	    else if (new_SNR > old_SNR)
	      {
		d_SNR_output[pos_y + pos_x] = new_SNR;
		d_SNR_taps[pos_y + pos_x] = nTaps;
	      }
	  }
	else
	  {
	    if (pos_x >= ( nColumns - offset ) && pos_x < nColumns)
	      {
		d_SNR_output[pos_y + pos_x] = 0;
		d_SNR_taps[pos_y + pos_x] = 0;
	      }
	  }
	pos_x = pos_x + WARP;
      }

  } //-------------------- KERNEL ENDS HERE --------------------------

  /** \brief Kernel wrapper function to SNR_GPU_limited kernel function. */
  void call_kernel_SNR_GPU_limited(const dim3 &grid_size, const dim3 &block_size, float *const d_FIR_input, float *const d_SNR_output,
				   ushort *const d_SNR_taps, float *const d_MSD, const int &x_steps, const int &nTaps,
				   const int &nColumns, const int &offset) {
    SNR_GPU_limited<<<grid_size, block_size>>>(d_FIR_input, d_SNR_output, d_SNR_taps, d_MSD,
					       x_steps, nTaps, nColumns, offset);
  }


} //namespace astroaccelerate
