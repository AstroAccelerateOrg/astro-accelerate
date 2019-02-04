#ifndef ASTRO_ACCELERATE_AA_DEVICE_SPS_INPLACE_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SPS_INPLACE_KERNEL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"

namespace astroaccelerate {

//This device constant is needed by the SPS (and SPS long) component and the SNR component
//It is also needed by device_load_data.cu and by device_single_pulse_search_kernel.cu
__device__ __constant__ float c_sqrt_taps[PD_MAXTAPS + 1];

void call_kernel_PD_ZC_GPU_KERNEL(const dim3 &grid_size, const dim3 &block_size, float *const d_input, float *const d_output,
				  const int &maxTaps, const int &nTimesamples, const int &nLoops);
void call_kernel_PD_INPLACE_GPU_KERNEL(const dim3 &grid_size, const dim3 &block_size, const int &SM_size, float *const d_input,
				       float *const d_temp, unsigned char *const d_output_taps, float *const d_MSD,
				       const int &maxTaps, const int &nTimesamples);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SPS_INPLACE_KERNEL_HPP
