// James Sharpe's peak finding code

#include "aa_device_peak_find_kernel.hpp"
#include "aa_device_peak_find_shared_kernel_functions.cuh"
#include "aa_device_threshold_shared_kernel_functions.cuh"
#include <stdio.h>

namespace astroaccelerate {

__global__ void peak_find_list2(const float *d_input, const int width, const int height, float threshold, int *gmem_pos, int shift, int DIT_value, ushort* d_input_taps, int max_peak_size, unsigned int *const d_peak_list_DM, unsigned int *const d_peak_list_TS, float *const d_peak_list_SNR, unsigned int *const d_peak_list_BW){

        int global_idx = blockDim.x*blockIdx.x + threadIdx.x;
        int list_pos;
	if (global_idx >= width*height) return;
        if (    (global_idx > 3*width - 1) && (global_idx < (height - 3)*width) &&
                ((global_idx % width > 2) && (global_idx % width < (width - 3)) ) ){
                float a[7];
                float b[7];
                float c[7];
                float d[7];
                float e[7];
                float f[7];
                float g[7];
                for (int i = 0; i < 7; i++){
                        a[i] = __ldg(d_input + global_idx - 3*width - 3 + i);
                        b[i] = __ldg(d_input + global_idx - 2*width - 3 + i);
                        c[i] = __ldg(d_input + global_idx - 1*width - 3 + i);
                        d[i] = __ldg(d_input + global_idx - 0*width - 3 + i);
                        e[i] = __ldg(d_input + global_idx + 1*width - 3 + i);
                        f[i] = __ldg(d_input + global_idx + 2*width - 3 + i);
                        g[i] = __ldg(d_input + global_idx + 3*width - 3 + i);
                }

                float center = d_input[global_idx];

                for (int i = 0; i < 7; i++){
                        if (center < a[i] ){
                                return;
                        } else if (center < b[i]){
                                return;
                        } else if (center < c[i]){
                                return;
                        } else if (center < e[i]){
                                return;
                        } else if (center < f[i]){
                                return;
                        } else if (center < g[i]){
                                return;
                        } else if (center < d[i]){
                                return;
			}
		}
                if (center > threshold){
                                        list_pos=atomicAdd(gmem_pos, 1);
                                        if(list_pos<max_peak_size){
						int h = (int)(global_idx/width);
						int wi = (global_idx % width);
//						printf("Center is : %lf %d %d %d %d %d\n", center, h, shift, height, global_idx, wi);
//						for (int i = 0; i < 7; i++){
//							printf("%4.1lf %4.1lf %4.1lf %4.1lf %4.1lf %4.1lf %4.1lf\n", a[i], b[i], c[i], d[i], e[i], f[i], g[i]);
//						}
                                                d_peak_list_DM[list_pos]  = h + shift; // DM coordinate (y)
                                                d_peak_list_TS[list_pos]  = wi*DIT_value + d_input_taps[global_idx]/2; // time coordinate (x)
                                                d_peak_list_SNR[list_pos] = center; // SNR value
                                                d_peak_list_BW[list_pos]  = d_input_taps[0]; // width of the boxcar
                                        } // if list_pos < max_peak
                                } // if threshold
        } // if outside boundaries
} // kernel end


  __global__ void dilate_peak_find(const float *d_input, ushort* d_input_taps,  unsigned int *d_peak_list_DM,  unsigned int *d_peak_list_TS, float *d_peak_list_SNR, unsigned int *d_peak_list_BW, const int width, const int height, const int offset, const float threshold, int max_peak_size, int *gmem_pos, int shift, int DIT_value) {
    int idxX = blockDim.x * blockIdx.x + threadIdx.x;
    int idxY = blockDim.y * blockIdx.y + threadIdx.y;
    if (idxX >= width-offset) return;
    if (idxY >= height) return;

    float dilated_value = 0.0f;
    float my_value = 0.0f;
    unsigned short peak = 0u;
    int list_pos;
    //handle boundary conditions - top edge
    if (idxY == 0) {
      //Special case for width of 1
      if (width == 1) {
	my_value = dilated_value = d_input[0];
      }
      //Top left corner case
      else if (idxX == 0) {
	float4 block = load_block_2x2(d_input, width);
	dilated_value = dilate4(block);
	my_value = block.x;
      } 
      //Top right corner case
      else if (idxX == (width-offset-1)) {
	float4 block = load_block_2x2(d_input+width-offset-2, width);
	dilated_value = dilate4(block);
	my_value = block.y;
      } else {
	float3x3 block = load_block_top(d_input, idxX, idxY, width);
	dilated_value = dilate3x3_top(block);
	my_value = block.y2;
      }
      //bottom edge
    } else if (idxY == height-1) {
      //Special case for width of 1
      if (width == 1) {
	my_value = dilated_value = d_input[width*(height-1)];
      }
      //Bottom left corner
      else if (idxX == 0) {
	float4 block = load_block_2x2(d_input+width*(height-2), width);
	dilated_value = dilate4(block);
	my_value = block.z;
      }
      //Bottom right corner
      else if (idxX == (width-offset-1)) {
	float4 block = load_block_2x2(d_input+width*(height-2)+width-offset-2, width);
	dilated_value = dilate4(block);
	my_value = block.w;
      } else {
	float3x3 block = load_block_bottom(d_input, idxX, idxY, width);        
	dilated_value = dilate3x3_bottom(block);
	my_value = block.y2;
      }
      //Left edge
    } else if (idxX == 0) {
      float3x3 block = load_block_left(d_input, idxX, idxY, width);        
      dilated_value = dilate3x3_left(block);
      my_value = block.y2;
    
      //right edge
    } else if (idxX == (width-offset-1)) {
      float3x3 block = load_block_right(d_input, idxX, idxY, width);        
      dilated_value = dilate3x3_right(block);
      my_value = block.y2;

    } else {
      float3x3 block = load_block(d_input, idxX, idxY, width);
      dilated_value = dilate3x3(block);
      my_value = block.y2;
    }
    peak = is_peak( my_value, dilated_value, threshold);
    if(peak==1){
      list_pos=atomicAdd(gmem_pos, 1);
      if(list_pos<max_peak_size){
	d_peak_list_DM[list_pos]  = idxY + shift; // DM coordinate (y)
	d_peak_list_TS[list_pos]  = idxX*DIT_value + d_input_taps[idxY*width+idxX]/2; // time coordinate (x)
	d_peak_list_SNR[list_pos] = my_value; // SNR value
	d_peak_list_BW[list_pos]  = d_input_taps[idxY*width+idxX]; // width of the boxcar
      }
    }
	
    //d_output[idxY*width+idxX] = peak;
  }


  __global__ void dilate_peak_find_for_fdas(const float *d_input, float *d_peak_list, float *d_MSD, const int width, const int height, const int offset, const float threshold, unsigned int max_peak_size, unsigned int *gmem_pos, float DM_trial) {
    int idxX = blockDim.x * blockIdx.x + threadIdx.x;
    int idxY = blockDim.y * blockIdx.y + threadIdx.y;
    if (idxX >= width-offset) return;
    if (idxY >= height) return;

    float dilated_value = 0.0f;
    float my_value = 0.0f;
    float SNR;
    int list_pos;
    //handle boundary conditions - top edge
    if (idxY == 0) {
      //Special case for width of 1
      if (width == 1) {
	my_value = dilated_value = d_input[0];
      }
      //Top left corner case
      else if (idxX == 0) {
	float4 block = load_block_2x2(d_input, width);
	dilated_value = dilate4(block);
	my_value = block.x;
      } 
      //Top right corner case
      else if (idxX == (width-offset-1)) {
	float4 block = load_block_2x2(d_input+width-offset-2, width);
	dilated_value = dilate4(block);
	my_value = block.y;
      } else {
	float3x3 block = load_block_top(d_input, idxX, idxY, width);
	dilated_value = dilate3x3_top(block);
	my_value = block.y2;
      }
      //bottom edge
    } else if (idxY == height-1) {
      //Special case for width of 1
      if (width == 1) {
	my_value = dilated_value = d_input[width*(height-1)];
      }
      //Bottom left corner
      else if (idxX == 0) {
	float4 block = load_block_2x2(d_input+width*(height-2), width);
	dilated_value = dilate4(block);
	my_value = block.z;
      }
      //Bottom right corner
      else if (idxX == (width-offset-1)) {
	float4 block = load_block_2x2(d_input+width*(height-2)+width-offset-2, width);
	dilated_value = dilate4(block);
	my_value = block.w;
      } else {
	float3x3 block = load_block_bottom(d_input, idxX, idxY, width);        
	dilated_value = dilate3x3_bottom(block);
	my_value = block.y2;
      }
      //Left edge
    } else if (idxX == 0) {
      float3x3 block = load_block_left(d_input, idxX, idxY, width);        
      dilated_value = dilate3x3_left(block);
      my_value = block.y2;
    
      //right edge
    } else if (idxX == (width-offset-1)) {
      float3x3 block = load_block_right(d_input, idxX, idxY, width);        
      dilated_value = dilate3x3_right(block);
      my_value = block.y2;

    } else {
      float3x3 block = load_block(d_input, idxX, idxY, width);
      dilated_value = dilate3x3(block);
      my_value = block.y2;
    }
    
    if(my_value == dilated_value){ // this means there is a peak
      SNR = (my_value-d_MSD[0])/d_MSD[1]; // calculation of SNR
      if(SNR > threshold) {
	list_pos=atomicAdd(gmem_pos, 1);
	if(list_pos<max_peak_size){
	  d_peak_list[4*list_pos]   = idxY;
	  d_peak_list[4*list_pos+1] = idxX;
	  d_peak_list[4*list_pos+2] = my_value;
	  d_peak_list[4*list_pos+3] = DM_trial;
	}
      }
    }
	
    //d_output[idxY*width+idxX] = peak;
  }


  // width DM
  // height time
  __global__ void dilate_peak_find_for_periods_old(const float *d_input, ushort* d_input_taps, float *d_peak_list, const int width, const int height, const int offset, const float threshold, int max_peak_size, int *gmem_pos, float *d_MSD, int DM_shift, int DIT_value) {
    int idxX = blockDim.x * blockIdx.x + threadIdx.x;
    int idxY = blockDim.y * blockIdx.y + threadIdx.y;
    if (idxX >= width-offset) return;
    if (idxY >= height) return;

    float dilated_value = 0.0f;
    float my_value = 0.0f;
    unsigned short peak = 0u;
    int list_pos;
    float mean = d_MSD[0];
    float sd   = d_MSD[1];
    //handle boundary conditions - top edge
    if (idxY == 0) {
      //Special case for width of 1
      if (width == 1) {
	my_value = dilated_value = d_input[0];
      }
      //Top left corner case
      else if (idxX == 0) {
	float4 block = load_block_2x2(d_input, width);
	dilated_value = dilate4(block);
	my_value = block.x;
      } 
      //Top right corner case
      else if (idxX == (width-offset-1)) {
	float4 block = load_block_2x2(d_input+width-offset-2, width);
	dilated_value = dilate4(block);
	my_value = block.y;
      } else {
	float3x3 block = load_block_top(d_input, idxX, idxY, width);
	dilated_value = dilate3x3_top(block);
	my_value = block.y2;
      }
      //bottom edge
    } else if (idxY == height-1) {
      //Special case for width of 1
      if (width == 1) {
	my_value = dilated_value = d_input[width*(height-1)];
      }
      //Bottom left corner
      else if (idxX == 0) {
	float4 block = load_block_2x2(d_input+width*(height-2), width);
	dilated_value = dilate4(block);
	my_value = block.z;
      }
      //Bottom right corner
      else if (idxX == (width-offset-1)) {
	float4 block = load_block_2x2(d_input+width*(height-2)+width-offset-2, width);
	dilated_value = dilate4(block);
	my_value = block.w;
      } else {
	float3x3 block = load_block_bottom(d_input, idxX, idxY, width);        
	dilated_value = dilate3x3_bottom(block);
	my_value = block.y2;
      }
      //Left edge
    } else if (idxX == 0) {
      float3x3 block = load_block_left(d_input, idxX, idxY, width);        
      dilated_value = dilate3x3_left(block);
      my_value = block.y2;
    
      //right edge
    } else if (idxX == (width-offset-1)) {
      float3x3 block = load_block_right(d_input, idxX, idxY, width);        
      dilated_value = dilate3x3_right(block);
      my_value = block.y2;

    } else {
      float3x3 block = load_block(d_input, idxX, idxY, width);
      dilated_value = dilate3x3(block);
      my_value = block.y2;
    }
    peak = is_peak( my_value, dilated_value, threshold);
    if(peak==1){
      list_pos=atomicAdd(gmem_pos, 1);
      if(list_pos<max_peak_size){
	d_peak_list[4*list_pos]   = idxX + DM_shift; // DM coordinate (y)
	d_peak_list[4*list_pos+1] = idxY/DIT_value; // frequency coordinate (x)
	float hrms = d_input_taps[idxY*width+idxX];
	d_peak_list[4*list_pos+3] = hrms; // width of the boxcar
	d_peak_list[4*list_pos+2] = inverse_white_noise(&my_value,&hrms,&mean,&sd); // SNR value
			
      }
    }
	
    //d_output[idxY*width+idxX] = peak;
  }






  __global__ void dilate_peak_find_for_periods(const float *d_input, ushort* d_input_taps, float *d_peak_list, const int width, const int height, const int offset, const float threshold, int max_peak_size, int *gmem_pos, float const* __restrict__ d_MSD, int DM_shift, int DIT_value) {
    int idxX = blockDim.x * blockIdx.x + threadIdx.x;
    int idxY = blockDim.y * blockIdx.y + threadIdx.y;
    if (idxX >= width-offset) return;
    if (idxY >= height) return;

    float dilated_value = 0.0f;
    float my_value = 0.0f;
    unsigned short peak = 0u;
    int list_pos;
    //handle boundary conditions - top edge
    if (idxY == 0) {
      //Special case for width of 1
      if (width == 1) {
	my_value = dilated_value = d_input[0];
      }
      //Top left corner case
      else if (idxX == 0) {
	float4 block = load_block_2x2(d_input, width);
	dilated_value = dilate4(block);
	my_value = block.x;
      } 
      //Top right corner case
      else if (idxX == (width-offset-1)) {
	float4 block = load_block_2x2(d_input+width-offset-2, width);
	dilated_value = dilate4(block);
	my_value = block.y;
      } else {
	float3x3 block = load_block_top(d_input, idxX, idxY, width);
	dilated_value = dilate3x3_top(block);
	my_value = block.y2;
      }
      //bottom edge
    } else if (idxY == height-1) {
      //Special case for width of 1
      if (width == 1) {
	my_value = dilated_value = d_input[width*(height-1)];
      }
      //Bottom left corner
      else if (idxX == 0) {
	float4 block = load_block_2x2(d_input+width*(height-2), width);
	dilated_value = dilate4(block);
	my_value = block.z;
      }
      //Bottom right corner
      else if (idxX == (width-offset-1)) {
	float4 block = load_block_2x2(d_input+width*(height-2)+width-offset-2, width);
	dilated_value = dilate4(block);
	my_value = block.w;
      } else {
	float3x3 block = load_block_bottom(d_input, idxX, idxY, width);        
	dilated_value = dilate3x3_bottom(block);
	my_value = block.y2;
      }
      //Left edge
    } else if (idxX == 0) {
      float3x3 block = load_block_left(d_input, idxX, idxY, width);        
      dilated_value = dilate3x3_left(block);
      my_value = block.y2;
    
      //right edge
    } else if (idxX == (width-offset-1)) {
      float3x3 block = load_block_right(d_input, idxX, idxY, width);        
      dilated_value = dilate3x3_right(block);
      my_value = block.y2;

    } else {
      float3x3 block = load_block(d_input, idxX, idxY, width);
      dilated_value = dilate3x3(block);
      my_value = block.y2;
    }
    peak = is_peak( my_value, dilated_value, threshold);
    if(peak==1){
      list_pos=atomicAdd(gmem_pos, 1);
      if(list_pos<max_peak_size){
	d_peak_list[4*list_pos]   = idxX + DM_shift; // DM coordinate (y)
	d_peak_list[4*list_pos+1] = idxY/DIT_value; // frequency coordinate (x)
	int hrms = (int) d_input_taps[idxY*width+idxX];
	d_peak_list[4*list_pos+3] = hrms; // width of the boxcar
	d_peak_list[4*list_pos+2] = my_value*__ldg(&d_MSD[2*hrms+1]) + __ldg(&d_MSD[2*hrms]); // SNR value
			
      }
    }
	
    //d_output[idxY*width+idxX] = peak;
  }


#define PPF_NTHREADS 64 
__global__ void gpu_Filter_peaks_kernel(unsigned int *d_new_peak_list_DM, unsigned int *d_new_peak_list_TS, unsigned int *d_new_peak_list_BW, float *d_new_peak_list_SNR, 
					unsigned int *d_peak_list_DM, unsigned int *d_peak_list_TS, unsigned int *d_peak_list_BW, float *d_peak_list_SNR,
					unsigned int nElements, unsigned int max_distance, int nLoops, int max_list_pos, int *gmem_pos){
	// PPF_DPB = 128 //this is because I set nThreads to 64
	// PPF_PEAKS_PER_BLOCK = something small like 10
	__shared__ float s_data_snr[PPF_DPB];
	__shared__ int s_data_dm[PPF_DPB];
	__shared__ int s_data_ts[PPF_DPB];
	__shared__ int s_flag[PPF_NTHREADS];
	int d, s;
	int elements_pos, pos;
	float snr, distance, fs, fd;
//	float4 f4temp;
	

	if(threadIdx.x<PPF_PEAKS_PER_BLOCK){
		s_flag[threadIdx.x] = 1;
	}
	else{
		s_flag[threadIdx.x] = 0;
	}

	
	for(int f=0; f<nLoops; f++){
		// Load new data blob
		//s_data[threadIdx.x + 2*PPF_DPB] = 0; // SNR
		//s_data[threadIdx.x + 64 + 2*PPF_DPB] = 0; // SNR
		
		pos = PPF_DPB*f + threadIdx.x;
		if(pos < nElements){
//			f4temp = __ldg(&d_peak_list[pos]);
			s_data_dm[threadIdx.x]  = d_peak_list_DM[pos]; //f4temp.x; // DM
			s_data_ts[threadIdx.x]  = d_peak_list_TS[pos]; //f4temp.y; // Time
			s_data_snr[threadIdx.x] = d_peak_list_SNR[pos]; //f4temp.z; // SNR
		}
	       	else {
                        s_data_dm[threadIdx.x]  = 0; //f4temp.x; // DM
                        s_data_ts[threadIdx.x]  = 0; //f4temp.y; // Time
                        s_data_snr[threadIdx.x] = -1000; //f4temp.z; // SNR
		}
//		if(blockIdx.x==0 && threadIdx.x==0) printf("point: [%d;%d;%lf]\n",  s_data_dm[threadIdx.x], s_data_ts[threadIdx.x], s_data_snr[threadIdx.x]);
		
		
		pos = PPF_DPB*f + threadIdx.x + PPF_NTHREADS;
		if(pos < nElements){
//			f4temp = __ldg(&d_peak_list[PPF_DPB*f + threadIdx.x + (PPF_DPB>>1)]);
			s_data_dm[threadIdx.x + PPF_NTHREADS ] = d_peak_list_DM[pos]; //f4temp.x; // DM
			s_data_ts[threadIdx.x + PPF_NTHREADS ] = d_peak_list_TS[pos]; //f4temp.y; // Time
			s_data_snr[threadIdx.x + PPF_NTHREADS] = d_peak_list_SNR[pos]; //f4temp.z; // SNR
		}
		else {
                        s_data_dm[threadIdx.x + PPF_NTHREADS]  = 0; //f4temp.x; // DM
                        s_data_ts[threadIdx.x + PPF_NTHREADS]  = 0; //f4temp.y; // Time
                        s_data_snr[threadIdx.x + PPF_NTHREADS] = -1000; //f4temp.z; // SNR		
		}
		
		__syncthreads();

		elements_pos = blockIdx.x*PPF_PEAKS_PER_BLOCK;
		for(int p=0; p<PPF_PEAKS_PER_BLOCK; p++){
//			if (blockIdx.x == 0) printf("%d %d\n", p, s_flag[p]);
			if((s_flag[p]) && ((elements_pos + p) < nElements)){
				//pos = elements_pos+p;
				//if(pos<nElements){
					d   = d_peak_list_DM[elements_pos+p]; // DM
					s   = d_peak_list_TS[elements_pos+p]; // Time
					snr = d_peak_list_SNR[elements_pos+p]; // SNR
					
					// first element
//					if(blockIdx.x==0) printf("s_data: %lf, snr: %lf, pos: %d\n", s_data_snr[threadIdx.x], snr, p);
					if( (s_data_snr[threadIdx.x] >= snr)){
						fs = ((float)s_data_dm[threadIdx.x] - (float)d);
						fd = ((float)s_data_ts[threadIdx.x] - (float)s);
						distance = (fd*fd + fs*fs);
//						if(blockIdx.x==0) printf("%d - %d = %d; %d - %d = %d\n",s_data_dm[threadIdx.x], d, fs, s_data_ts[threadIdx.x], s, fd, distance);
						if( (distance < (float)max_distance) && (distance!=0) ){
//							if(blockIdx.x==0) printf("distance: %d %lf %lf %lf %d %d;\n", p, distance, fs, fd, s, d);
							s_flag[p]=0;
						}
					}
					
					//second element
					if(s_data_snr[threadIdx.x + PPF_NTHREADS] >= snr){
						fs = ((float)s_data_dm[threadIdx.x + PPF_NTHREADS] - (float)d);
						fd = ((float)s_data_ts[threadIdx.x + PPF_NTHREADS] - (float)s);
						distance = (fd*fd + fs*fs);
//						if(blockIdx.x==0) printf("%d - %d = %d; %d - %d = %d\n",s_data_dm[threadIdx.x], d, fs, s_data_ts[threadIdx.x], s, fd, distance);
						if( (distance < (float)max_distance) && (distance!=0)){
							s_flag[p]=0;
//							if(blockIdx.x==0) printf("xdistance: %d %lf %lf %lf %d %d;\n", p, distance, fs, fd, s, d);
						}
					}
				//}
			}
		} // for p
		
	}
	
	// Saving peaks that got through
	elements_pos = blockIdx.x*PPF_PEAKS_PER_BLOCK;
	if(threadIdx.x < PPF_PEAKS_PER_BLOCK){
		if( (s_flag[threadIdx.x] == 1) && ((elements_pos + threadIdx.x) < nElements)){
			int list_pos=atomicAdd(gmem_pos, 1);
			if(list_pos<max_list_pos){
				d_new_peak_list_DM[list_pos]  = d_peak_list_DM[elements_pos  + threadIdx.x];
				d_new_peak_list_TS[list_pos]  = d_peak_list_TS[elements_pos  + threadIdx.x];
				d_new_peak_list_BW[list_pos]  = d_peak_list_BW[elements_pos  + threadIdx.x];
				d_new_peak_list_SNR[list_pos] = d_peak_list_SNR[elements_pos + threadIdx.x];
			}
		}
	}
}







  /** \brief Kernel wrapper function for dilate_peak_find kernel function. */
  void call_kernel_dilate_peak_find(const dim3 &grid_size, const dim3 &block_size,
				    float *const d_input, ushort *const d_input_taps,  unsigned int *const d_peak_list_DM,
				    unsigned int *const d_peak_list_TS, float *const d_peak_list_SNR, unsigned int *const d_peak_list_BW,
				    const int &width, const int &height, const int &offset, const float &threshold,
				    const int &max_peak_size, int *const gmem_pos, const int &shift, const int &DIT_value) {
    dilate_peak_find<<<grid_size, block_size>>>(d_input, d_input_taps, d_peak_list_DM,
						d_peak_list_TS, d_peak_list_SNR, d_peak_list_BW,
						width, height, offset, threshold,
						max_peak_size, gmem_pos, shift, DIT_value);
  }

  /** \brief Kernel wrapper function for dilate_peak_find_for_fdas kernel function. */
  void call_kernel_dilate_peak_find_for_fdas(const dim3 &grid_size, const dim3 &block_size,
					     float *const d_input, float *const d_peak_list, float *const d_MSD, const int &width,
					     const int &height, const int &offset, const float &threshold,
					     const unsigned int &max_peak_size, unsigned int *const gmem_pos, const float &DM_trial) {
    dilate_peak_find_for_fdas<<<grid_size, block_size>>>( d_input, d_peak_list, d_MSD, width,
							  height, offset, threshold,
							  max_peak_size, gmem_pos, DM_trial);
  }

  /** \brief Kernel wrapper function for dilate_peak_find_for_periods_old kernel function. */
  void call_kernel_dilate_peak_find_for_periods_old(const dim3 &grid_size, const dim3 &block_size,
						    float *const d_input, ushort *const d_input_taps, float *const d_peak_list,
						    const int &width, const int &height, const int &offset,
						    const float &threshold,
						    const int &max_peak_size, int *const gmem_pos, float *const d_MDF,
						    const int &DM_shift, const int &DIT_value) {
    dilate_peak_find_for_periods_old<<<grid_size, block_size>>>(d_input, d_input_taps, d_peak_list,
								width, height, offset,
								threshold,
								max_peak_size, gmem_pos, d_MDF,
								DM_shift, DIT_value);
  }

  /** \brief Kernel wrapper function for _dilate_peak_find_for_periods kernel function. */
  void call_kernel_dilate_peak_find_for_periods(const dim3 &grid_size, const dim3 &block_size,
						float *const d_input, ushort *const d_input_taps,
						float *const d_peak_list, const int &width,
						const int &height, const int &offset, const float &threshold,
						const int &max_peak_size,
						int *const gmem_pos, float const *const d_MSD, const int &DM_shift, const int &DIT_value) {
    dilate_peak_find_for_periods<<<grid_size, block_size>>>(d_input, d_input_taps,
							    d_peak_list, width,
							    height, offset, threshold,
							    max_peak_size,
							    gmem_pos, d_MSD, DM_shift, DIT_value);
  }

	void call_kernel_peak_find_list(const dim3 &grid_size, const dim3 &block_size, float *const d_input, const int width, const int height, const float &threshold, int *const gmem_pos, const int &shift, const int &DIT_value, ushort *const d_input_taps, const int &max_peak_size, unsigned int *const d_peak_list_DM, unsigned int *const d_peak_list_TS, float *const d_peak_list_SNR, unsigned int *const d_peak_list_BW){
		peak_find_list2<<<grid_size,block_size>>>(d_input, width, height, threshold, gmem_pos, shift, DIT_value, d_input_taps, max_peak_size, d_peak_list_DM, d_peak_list_TS, d_peak_list_SNR, d_peak_list_BW);
	}

	void call_gpu_Filter_peaks(unsigned int *new_peak_list_DM, unsigned int *new_peak_list_TS, unsigned int *new_peak_list_BW, float *new_peak_list_SNR, unsigned int *d_peak_list_DM, unsigned int *d_peak_list_TS, unsigned int *d_peak_list_BW, float *d_peak_list_SNR, unsigned int nElements, unsigned int max_distance, int max_list_pos, int *gmem_pos){
		int nThreads, nBlocks_x, nLoops;
		nThreads = 64;
		nBlocks_x = (nElements + PPF_PEAKS_PER_BLOCK - 1)/PPF_PEAKS_PER_BLOCK; //10
		nLoops = (nElements + PPF_DPB - 1)/PPF_DPB; //128

		dim3 blockDim(nThreads, 1, 1);
		dim3 gridSize(nBlocks_x, 1, 1);
		if ( (nBlocks_x > 0) && (nLoops > 0) ){
			gpu_Filter_peaks_kernel<<<gridSize, blockDim>>>(new_peak_list_DM, new_peak_list_TS, new_peak_list_BW, new_peak_list_SNR, d_peak_list_DM, d_peak_list_TS, d_peak_list_BW, d_peak_list_SNR, nElements, max_distance*max_distance, nLoops, max_list_pos, gmem_pos);
		}
	}

} //namespace astroaccelerate
