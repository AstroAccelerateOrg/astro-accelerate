#ifndef ASTRO_ACCELERATE_DEVICE_THRESHOLD_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_THRESHOLD_KERNEL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include "params.hpp"

void call_kernel_THR_GPU_WARP(dim3 grid_size, dim3 block_size,
			      float const* d_input, ushort *d_input_taps, unsigned int *d_output_list_DM,
			      unsigned int *d_output_list_TS, float *d_output_list_SNR,
			      unsigned int *d_output_list_BW, int *gmem_pos,
			      float threshold, int nTimesamples, int offset, int shift, int max_list_size,
			      int DIT_value);

void call_kernel_GPU_Threshold_for_periodicity_kernel_old(dim3 grid_size, dim3 block_size,
							  float const* d_input, ushort *d_input_harms,
							  float *d_output_list, int *gmem_pos, float *d_MSD,
							  float threshold, int primary_size,
							  int secondary_size, int DM_shift,
							  int max_list_size, int DIT_value);

void call_kernel_GPU_Threshold_for_periodicity_kernel(dim3 grid_size, dim3 block_size,
						      float const*  d_input, ushort *d_input_harms,
						      float *d_output_list, int *gmem_pos, float const* d_MSD,
						      float threshold, int primary_size, int secondary_size,
						      int DM_shift, int max_list_size, int DIT_value);



#endif
