// Added by Karel Adamek

#include "aa_device_threshold_kernel.hpp"
#include "aa_device_peak_find_shared_kernel_functions.cuh"
#include "aa_device_threshold_shared_kernel_functions.cuh"
#include "aa_device_cuda_deprecated_wrappers.cuh"

namespace astroaccelerate {

  __global__ void THR_GPU_WARP(float const* __restrict__ d_input, ushort *d_input_taps, unsigned int *d_output_list_DM, unsigned int *d_output_list_TS, float *d_output_list_SNR, unsigned int *d_output_list_BW, unsigned int *gmem_pos, float threshold, int nTimesamples, int offset, int shift, unsigned int max_list_size, int DIT_value) {
    int local_id;
    local_id = threadIdx.x & (WARP - 1);
    int warp_id;
    warp_id = threadIdx.x>>5;
	
    int pos_x, pos_y, mask, leader;
    unsigned int list_pos;
    float R;
	
    pos_y = (blockIdx.y*THR_WARPS_PER_BLOCK + warp_id)*nTimesamples;
    pos_x = blockIdx.x*WARP*THR_ELEM_PER_THREAD + local_id;
	
    for(int i=0; i<THR_ELEM_PER_THREAD; i++){
      if(pos_x<offset){
	R=__ldg(&d_input[pos_x + pos_y]);
	if(R > threshold) {
	  mask=aa_ballot(AA_ASSUME_MASK,1);
	  leader=__ffs(mask)-1;
	  if(local_id==leader) list_pos=atomicAdd(gmem_pos,__popc(mask));
	  list_pos=aa_shfl(AA_ASSUME_MASK,list_pos,leader);
	  list_pos=list_pos+__popc(mask&((1<<local_id)-1));
	  if(list_pos<max_list_size){
	    d_output_list_DM[list_pos]  = blockIdx.y*THR_WARPS_PER_BLOCK + warp_id + shift;
	    d_output_list_TS[list_pos]  = DIT_value*pos_x + d_input_taps[(pos_x + pos_y)]/2;
	    d_output_list_SNR[list_pos] = R;
	    d_output_list_BW[list_pos]  = (unsigned int) d_input_taps[(pos_x + pos_y)];
	  }
	}
      }
      pos_x = pos_x + WARP;
    }
  }

  __global__ void GPU_Threshold_for_periodicity_normal_kernel(float const* __restrict__ d_input, ushort *d_input_harms, float *d_output_list, unsigned int *gmem_pos, float const* __restrict__ d_MSD, float threshold, int nTimesamples, int nDMs, int DM_shift, unsigned int max_list_size, int DIT_value) {
    int pos_x, pos_y, pos, mask, leader;
    unsigned int list_pos;
    float R;
    int hrms;
    
    pos_x = blockIdx.x*blockDim.x + threadIdx.x; // timesamples
    pos_y = blockIdx.y*THR_ELEM_PER_THREAD; // DM-trials
    
    for(int i=0; i<THR_ELEM_PER_THREAD; i++){
        if( (pos_x<nTimesamples) && (pos_y<nDMs)){
            pos = pos_y*nTimesamples + pos_x;
            
            //--------> Thresholding
            R = __ldg(&d_input[pos]);
            if(R > threshold) {
                mask=aa_ballot(AA_ASSUME_MASK,1);
                leader=__ffs(mask)-1;
                if(threadIdx.x==leader) list_pos=atomicAdd(gmem_pos,__popc(mask));
                list_pos=aa_shfl(AA_ASSUME_MASK,list_pos,leader);
                list_pos=list_pos+__popc(mask&((1<<threadIdx.x)-1));
                if(list_pos<max_list_size){
                    d_output_list[4*list_pos]   = pos_y + DM_shift;
                    d_output_list[4*list_pos+1] = pos_x/DIT_value;
                    hrms = (int) d_input_harms[pos];
                    d_output_list[4*list_pos+3] = (float) hrms;
                    d_output_list[4*list_pos+2] = R*__ldg(&d_MSD[2*hrms+1]) + __ldg(&d_MSD[2*hrms]);
                }
            }
            //-------------------------<
        }
        pos_y++;
    }
  }

  __global__ void GPU_Threshold_for_periodicity_transposed_kernel(float const* __restrict__ d_input, ushort *d_input_harms, float *d_output_list, unsigned int *gmem_pos, float const* __restrict__ d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, unsigned int max_list_size, int DIT_value) {
    int pos_p, pos_s, pos, mask, leader;
    unsigned int list_pos;
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
          mask=aa_ballot(AA_ASSUME_MASK,1);
          leader=__ffs(mask)-1;
          if(threadIdx.x==leader) list_pos=atomicAdd(gmem_pos,__popc(mask));
          list_pos=aa_shfl(AA_ASSUME_MASK,list_pos,leader);
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

  /** \brief Kernel wrapper function for THR_GPU_WARP kernel function. */
  void call_kernel_THR_GPU_WARP(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      ushort *const d_input_taps,
      unsigned int *const d_output_list_DM,
      unsigned int *const d_output_list_TS,
      float *const d_output_list_SNR,
      unsigned int *const d_output_list_BW,
      unsigned int *const gmem_pos,
      const float &threshold,
      const int &nTimesamples,
      const int &offset,
      const int &shift,
      const unsigned int &max_list_size,
      const int &DIT_value
  ) {
      THR_GPU_WARP<<<grid_size, block_size>>>(
          d_input,
          d_input_taps,
          d_output_list_DM,
          d_output_list_TS,
          d_output_list_SNR,
          d_output_list_BW,
          gmem_pos,
          threshold,
          nTimesamples,
          offset,
          shift,
          max_list_size,
          DIT_value
      );
  }

  /** \brief Kernel wrapper function for GPU_Threshold_for_periodicity_kernel_old kernel function. */
  void call_kernel_GPU_Threshold_for_periodicity_normal_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      ushort *const d_input_harms,
      float *const d_output_list,
      unsigned int *const gmem_pos,
      float const *const d_MSD,
      const float &threshold,
      const int &nTimesamples,
      const int &nDMs,
      const int &DM_shift,
      const unsigned int &max_list_size,
      const int &DIT_value
  ) {
      GPU_Threshold_for_periodicity_normal_kernel<<<grid_size, block_size>>>(
          d_input,
          d_input_harms,
          d_output_list,
          gmem_pos,
          d_MSD,
          threshold,
          nTimesamples,
          nDMs,
          DM_shift,
          max_list_size,
          DIT_value
      );

  }

  /** \brief Kernel wrapper function for GPU_Threshold_for_periodicity_kernel kernel function. */
  void call_kernel_GPU_Threshold_for_periodicity_transposed_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      ushort *const d_input_harms,
      float *const d_output_list,
      unsigned int *const gmem_pos,
      float const *const d_MSD,
      const float &threshold,
      const int &primary_size,
      const int &secondary_size,
      const int &DM_shift,
      const unsigned int &max_list_size,
      const int &DIT_value
  ) {
      GPU_Threshold_for_periodicity_transposed_kernel<<<grid_size, block_size>>>(
          d_input,
          d_input_harms,
          d_output_list,
          gmem_pos,
          d_MSD,
          threshold,
          primary_size,
          secondary_size,
          DM_shift,
          max_list_size,
          DIT_value
      );
  }
} //namespace astroaccelerate
