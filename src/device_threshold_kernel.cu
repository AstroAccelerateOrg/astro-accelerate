// Added by Karel Adamek

#include "device_threshold_kernel.hpp"
#include "device_peak_find_shared_kernel_functions.cuh"
#include "device_threshold_shared_kernel_functions.cuh"

__global__ void THR_GPU_WARP(float const* __restrict__ d_input, ushort *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nTimesamples, int offset, int shift, int max_list_size, int DIT_value, float dm_step, float dm_low, float sampling_time, float inBin, float start_time) {
	int local_id;
	local_id = threadIdx.x & (WARP - 1);
	int warp_id;
	warp_id = threadIdx.x>>5;
	
	int pos_x, pos_y, list_pos, mask, leader;
	float R;

	pos_y = (blockIdx.y*THR_WARPS_PER_BLOCK + warp_id)*nTimesamples;
	pos_x = blockIdx.x*WARP*THR_ELEM_PER_THREAD + local_id;
	
	for(int i=0; i<THR_ELEM_PER_THREAD; i++){
		if(pos_x<offset){
			R=__ldg(&d_input[pos_x + pos_y]);
			if(R > threshold) {
				mask=__ballot(1);
				leader=__ffs(mask)-1;
				if(local_id==leader) list_pos=atomicAdd(gmem_pos,__popc(mask));
				list_pos=__shfl(list_pos,leader);
				list_pos=list_pos+__popc(mask&((1<<local_id)-1));
				if(list_pos<max_list_size){
					d_output_list[4*list_pos]   = (blockIdx.y*THR_WARPS_PER_BLOCK + warp_id + shift)*dm_step + dm_low;
					d_output_list[4*list_pos+1] = (DIT_value*pos_x + (float) d_input_taps[(pos_x + pos_y)]/2.0)*sampling_time + start_time;
					d_output_list[4*list_pos+2] = R;
					d_output_list[4*list_pos+3] = ((float) d_input_taps[(pos_x + pos_y)])*inBin;
				}
			}
		}
		pos_x = pos_x + WARP;
	}
}

__global__ void GPU_Threshold_for_periodicity_kernel_old(float const* __restrict__ d_input, ushort *d_input_harms, float *d_output_list, int *gmem_pos, float *d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int max_list_size, int DIT_value) {
	int pos_p, pos_s, pos, list_pos, mask, leader;
	float R;
	float hrms;
	float mean = d_MSD[0];
	float sd   = d_MSD[1];
	
	pos_p = blockIdx.x*blockDim.x*THR_ELEM_PER_THREAD + threadIdx.x;
	pos_s = blockIdx.y*blockDim.y + threadIdx.y;
	
	
	for(int i=0; i<THR_ELEM_PER_THREAD; i++){
		if( (pos_p<primary_size) && (pos_s<secondary_size)){
			pos = pos_s*primary_size + pos_p;
			
			//--------> Thresholding
			R = __ldg(&d_input[pos]);
			if(R > threshold) {
				mask=__ballot(1);
				leader=__ffs(mask)-1;
				if(threadIdx.x==leader) list_pos=atomicAdd(gmem_pos,__popc(mask));
				list_pos=__shfl(list_pos,leader);
				list_pos=list_pos+__popc(mask&((1<<threadIdx.x)-1));
				if(list_pos<max_list_size){
					d_output_list[4*list_pos]   = pos_p + DM_shift;
					d_output_list[4*list_pos+1] = pos_s/DIT_value;
					hrms = (float) d_input_harms[pos];
					d_output_list[4*list_pos+3] = hrms;
					d_output_list[4*list_pos+2] = inverse_white_noise(&R,&hrms,&mean,&sd);
				}
			}
			//-------------------------<
			
		}
		pos_p = pos_p + blockDim.x;
	}
}

__global__ void GPU_Threshold_for_periodicity_kernel(float const* __restrict__ d_input, ushort *d_input_harms, float *d_output_list, int *gmem_pos, float const* __restrict__ d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int max_list_size, int DIT_value) {
	int pos_p, pos_s, pos, list_pos, mask, leader;
	float R;
	int hrms;
	
	pos_p = blockIdx.x*blockDim.x*THR_ELEM_PER_THREAD + threadIdx.x;
	pos_s = blockIdx.y*blockDim.y + threadIdx.y;
	
	
	for(int i=0; i<THR_ELEM_PER_THREAD; i++){
		if( (pos_p<primary_size) && (pos_s<secondary_size)){
			pos = pos_s*primary_size + pos_p;
			
			//--------> Thresholding
			R = __ldg(&d_input[pos]);
			if(R > threshold) {
				mask=__ballot(1);
				leader=__ffs(mask)-1;
				if(threadIdx.x==leader) list_pos=atomicAdd(gmem_pos,__popc(mask));
				list_pos=__shfl(list_pos,leader);
				list_pos=list_pos+__popc(mask&((1<<threadIdx.x)-1));
				if(list_pos<max_list_size){
					d_output_list[4*list_pos]   = pos_p + DM_shift;
					d_output_list[4*list_pos+1] = pos_s/DIT_value;
					hrms = (int) d_input_harms[pos];
					d_output_list[4*list_pos+3] = hrms;
					d_output_list[4*list_pos+2] = R*__ldg(&d_MSD[2*hrms+1]) + __ldg(&d_MSD[2*hrms]);
				}
			}
			//-------------------------<
			
		}
		pos_p = pos_p + blockDim.x;
	}
}

void call_kernel_THR_GPU_WARP(dim3 grid_size, dim3 block_size,
			      float const* d_input, ushort *d_input_taps,
			      unsigned int *d_output_list_DM, unsigned int *d_output_list_TS,
			      float *d_output_list_SNR, unsigned int *d_output_list_BW,
			      int *gmem_pos, float threshold, int nTimesamples, int offset,
			      int shift, int max_list_size, int DIT_value) {
  THR_GPU_WARP<<<grid_size, block_size>>>(d_input, d_input_taps, d_output_list_DM, d_output_list_TS,
					  d_output_list_SNR, d_output_list_BW,
					  gmem_pos, threshold, nTimesamples, offset,
					  shift, max_list_size, DIT_value);
}

void call_kernel_GPU_Threshold_for_periodicity_kernel_old(dim3 grid_size, dim3 block_size,
							  float const* d_input, ushort *d_input_harms,
							  float *d_output_list, int *gmem_pos, float *d_MSD,
							  float threshold, int primary_size, int secondary_size,
							  int DM_shift, int max_list_size, int DIT_value) {
  GPU_Threshold_for_periodicity_kernel_old<<<grid_size, block_size>>>(d_input, d_input_harms, d_output_list,
								      gmem_pos, d_MSD, threshold, primary_size,
								      secondary_size, DM_shift, max_list_size, DIT_value);

}

void call_kernel_GPU_Threshold_for_periodicity_kernel(dim3 grid_size, dim3 block_size,
						      float const* d_input, ushort *d_input_harms,
						      float *d_output_list, int *gmem_pos, float const* d_MSD,
						      float threshold, int primary_size,
						      int secondary_size, int DM_shift, int max_list_size, int DIT_value) {
  GPU_Threshold_for_periodicity_kernel<<<grid_size, block_size>>>(d_input, d_input_harms,
								 d_output_list, gmem_pos, d_MSD,
								 threshold, primary_size,
								 secondary_size, DM_shift, max_list_size, DIT_value);
}
