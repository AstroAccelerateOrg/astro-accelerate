// James Sharpe's peak finding code

#ifndef PEAK_FIND_KERNEL_
#define PEAK_FIND_KERNEL_

struct float3x3 {
    float x1, x2, x3;
    float y1, y2, y3;
    float z1, z2, z3;
};

__device__ ushort is_peak(const float original_value, const float dilated_value, const float threshold) {
    return (original_value > threshold && original_value == dilated_value) ? 1u: 0u;
}

__device__ float3 load_row(const float *input, int idx) {
    float3 val;
    val.x = __ldg(input+idx);
    val.y = __ldg(input+idx+1);
    val.z = __ldg(input+idx+2);
    return val;
}

__device__ float2 load_row2(const float *input, int idx) {
    float2 val;
	val.x=__ldg(input+idx);
	val.y=__ldg(input+idx+1);
    return val;
}

__device__ float3x3 load_block(const float *input, int idxX, int idxY, int width) {
    int idx = idxY*width + idxX;
    float3x3 data;
    float3 yrow = load_row(input, idx-1);
    float3 xrow = load_row(input, idx-width-1);
    float3 zrow = load_row(input, idx+width-1);
    data.x1 = xrow.x;
    data.x2 = xrow.y;
    data.x3 = xrow.z;
    data.y1 = yrow.x;
    data.y2 = yrow.y;
    data.y3 = yrow.z;
    data.z1 = zrow.x;
    data.z2 = zrow.y;
    data.z3 = zrow.z;
    return data;
}

__device__ float3x3 load_block_top(const float *input, int idxX, int idxY, int width) {
     int idx = idxY*width + idxX;
     float3x3 data;
     float3 yrow = load_row(input, idx-1);
     float3 zrow = load_row(input, idx+width-1);
     data.y1 = yrow.x;
     data.y2 = yrow.y;
     data.y3 = yrow.z;
     data.z1 = zrow.x;
     data.z2 = zrow.y;
     data.z3 = zrow.z;
     return data;
}

__device__ float3x3 load_block_bottom(const float *input, int idxX, int idxY, int width) {
     int idx = idxY*width + idxX;
     float3x3 data;
     float3 yrow = load_row(input, idx-1);
     float3 xrow = load_row(input, idx-width-1);
     data.y1 = yrow.x;
     data.y2 = yrow.y;
     data.y3 = yrow.z;
     data.x1 = xrow.x;
     data.x2 = xrow.y;
     data.x3 = xrow.z;
     return data;
}

__device__ float3x3 load_block_left(const float *input, int idxX, int idxY, int width) {
     int idx = idxY*width + idxX;
     float3x3 data;
     float2 xrow = load_row2(input, idx-width);
     float2 yrow = load_row2(input, idx);
     float2 zrow = load_row2(input, idx+width);
     data.x2 = xrow.x;
     data.x3 = xrow.y;
     data.y2 = yrow.x;
     data.y3 = yrow.y;
     data.z2 = zrow.x;
     data.z3 = zrow.y;
     return data;
}

__device__ float3x3 load_block_right(const float *input, int idxX, int idxY, int width) {
     int idx = idxY*width + idxX;
     float3x3 data;
     float2 xrow = load_row2(input, idx-width-1);
     float2 yrow = load_row2(input, idx-1);
     float2 zrow = load_row2(input, idx+width-1);
     data.x1 = xrow.x;
     data.x2 = xrow.y;
     data.y1 = yrow.x;
     data.y2 = yrow.y;
     data.z1 = zrow.x;
     data.z2 = zrow.y;
     return data;
}

__device__ float4 load_block_2x2(const float *input, int width) {
	float2 first, second;
	first.x  = __ldg(input);        first.y  = __ldg(input+1);
	second.x = __ldg(input+width);  second.y = __ldg(input+width+1);
	float4 result;
	result.x = first.x; result.y = first.y; result.z = second.x; result.w = second.y;
	return result;
}



__device__ float dilate3x3_left(const float3x3 i) {
     float max = fmaxf(i.x2, i.y2);
     max = fmaxf(max, i.x3);
     max = fmaxf(max, i.y3);
     max = fmaxf(max, i.z2);
     max = fmaxf(max, i.z3);
     return max;
}

__device__ float dilate3x3_right(const float3x3 i) {
     float max = fmaxf(i.x2, i.y2);
     max = fmaxf(max, i.x1);
     max = fmaxf(max, i.y1);
     max = fmaxf(max, i.z2);
     max = fmaxf(max, i.z1);
     return max;
}

__device__ float dilate3x3_top(const float3x3 i) {
     float max = fmaxf(i.y1, i.y2);
     max = fmaxf(max, i.y3);
     max = fmaxf(max, i.z1);
     max = fmaxf(max, i.z2);
     max = fmaxf(max, i.z3);
     return max;
}

__device__ float dilate3x3_bottom(const float3x3 i) {
     float max = fmaxf(i.y1, i.y2);
     max = fmaxf(max, i.y3);
     max = fmaxf(max, i.x1);
     max = fmaxf(max, i.x2);
     max = fmaxf(max, i.x3);
     return max;
}

__device__ float dilate4(const float4 i) {
    float max = fmaxf(i.x, i.y);
    max = fmaxf(max, i.z);
    max = fmaxf(max, i.w);
    return max;
}

__device__ float dilate3x3(const float3x3 i) {
     float max = fmaxf(i.x1, i.x2);
     max = fmaxf(max, i.x3);
     max = fmaxf(max, i.y1);
     max = fmaxf(max, i.y2);
     max = fmaxf(max, i.y3);
     max = fmaxf(max, i.z1);
     max = fmaxf(max, i.z2);
     max = fmaxf(max, i.z3);
     return max;
}


__global__ void dilate_peak_find(const float *d_input, ushort* d_input_taps, float *d_peak_list, const int width, const int height, const int offset, const float threshold, int max_peak_size, int *gmem_pos, int shift, int DIT_value, float dm_step, float dm_low, float sampling_time, float inBin, float start_time) {
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
			d_peak_list[4*list_pos]   = (idxY + shift)*dm_step + dm_low; // DM coordinate (y)
			d_peak_list[4*list_pos+1] = (idxX*DIT_value + (float) d_input_taps[idxY*width+idxX]/2.0)*sampling_time + start_time; // time coordinate (x)
			d_peak_list[4*list_pos+2] = my_value; // SNE value
			d_peak_list[4*list_pos+3] = ((float) d_input_taps[idxY*width+idxX])*inBin; // width of the boxcar
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







#endif