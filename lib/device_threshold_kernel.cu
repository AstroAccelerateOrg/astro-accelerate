// Added by Karel Adamek 

#ifndef THRESHOLD_KERNEL_H_
#define THRESHOLD_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

__global__ void THR_GPU_WARP(float2 const* __restrict__ d_input, ushort *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nTimesamples, int offset, int shift, int max_list_size, int DIT_value) {
	int local_id;
	local_id = threadIdx.x & (WARP - 1);
	int warp_id;
	warp_id = threadIdx.x>>5;
	offset=offset>>1;
	
	int pos_x, pos_y, list_pos, mask, leader;
	float2 R;
	
	pos_y = (blockIdx.y*THR_WARPS_PER_BLOCK + warp_id)*(nTimesamples>>1);
	pos_x = blockIdx.x*WARP*THR_ELEM_PER_THREAD + local_id;
	
	for(int i=0; i<THR_ELEM_PER_THREAD; i++){
		if(pos_x<offset){
			R=__ldg(&d_input[pos_x + pos_y]);
			if(R.x > threshold) {
				mask=__ballot(1);
				leader=__ffs(mask)-1;
				if(local_id==leader) list_pos=atomicAdd(gmem_pos,__popc(mask));
				list_pos=__shfl(list_pos,leader);
				list_pos=list_pos+__popc(mask&((1<<local_id)-1));
				if(list_pos<max_list_size){
					d_output_list[4*list_pos]   = blockIdx.y*THR_WARPS_PER_BLOCK + warp_id + shift;
					d_output_list[4*list_pos+1] = 2*DIT_value*((blockIdx.x*THR_ELEM_PER_THREAD + i)*WARP + local_id);
					d_output_list[4*list_pos+2] = R.x;
					d_output_list[4*list_pos+3] = (float) d_input_taps[2*(pos_x + pos_y)];
				}
			}
			if(R.y > threshold) {
				mask=__ballot(1);
				leader=__ffs(mask)-1;
				if(local_id==leader) list_pos=atomicAdd(gmem_pos,__popc(mask));
				list_pos=__shfl(list_pos,leader);
				list_pos=list_pos+__popc(mask&((1<<local_id)-1));
				if(list_pos<max_list_size){
					d_output_list[4*list_pos]   = blockIdx.y*THR_WARPS_PER_BLOCK + warp_id + shift;
					d_output_list[4*list_pos+1] = 2*DIT_value*((blockIdx.x*THR_ELEM_PER_THREAD + i)*WARP + local_id)+1;
					d_output_list[4*list_pos+2] = R.y;
					d_output_list[4*list_pos+3] = (float) d_input_taps[2*(pos_x + pos_y) + 1];
				}
			}
		}
		pos_x = pos_x + WARP;
	}
}
#endif
