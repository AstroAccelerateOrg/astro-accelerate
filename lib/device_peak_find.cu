//-*-c++-*-

#include <npp.h>
#include "helper_cuda.h"

#include <stdint.h>
#include <cassert>
#include <vector>
#include <memory>
#include <iostream>

class NppImage {
public:
     NppImage(int32_t width, int32_t height, Npp32f *data=nullptr, int32_t linestep=0)
	: data(data), linestep(linestep), width(width), height(height)
     {
     }

     
     NppImage(const NppImage &) = delete;
     NppImage & operator = (const NppImage &) = delete;
     NppiSize size() const {
	return NppiSize() = {width, height};
     }

     Npp32f * data;
     Npp32s linestep; 
protected:
     /**
      * Only allow derived classes who can control the lifetime
      * of the memory handled in the data pointer to move
      */
     NppImage(NppImage && other) {
         this->data = other.data;
         this->width = other.width;
         this->height = other.height;
         this->linestep = other.linestep;
         other.data = nullptr;
      }
     Npp32s width, height;
};

class AllocatedNppImage : public NppImage {
public:
     AllocatedNppImage(int32_t width, int32_t height)
        : NppImage(width, height)
     {
         this->data = nppiMalloc_32f_C1(width, height, &this->linestep);
	 if (data == nullptr) throw;
     }
     
     AllocatedNppImage(AllocatedNppImage && other) 
	: NppImage(std::move(other))
     {
     }

     AllocatedNppImage & operator = (AllocatedNppImage && other) {
         this->data = other.data;
         this->width = other.width;
         this->height = other.height;
         this->linestep = other.linestep;
         other.data = nullptr;
         return *this;
     }

     ~AllocatedNppImage() {
         nppiFree(data);
     }
};

void handle_npp_error(NppStatus status)
{
    printf("NPP Error: %s", _cudaGetErrorEnum(status));
    exit(-1);
}

/**
 * Represents a 3x3 block of data laid out as follows:
 * 
 *  x1  x2  x3
 *  y1  y2  y3
 *  z1  z2  z3
 *
 * On the edges of the 2d array the values populated are
 * defined by the border handling of the load code
 */
struct float3x3 {
    float x1, x2, x3;
    float y1, y2, y3;
    float z1, z2, z3;
};

__device__ Npp16u is_peak(const float original_value, const float dilated_value, const float threshold=0.0f)
{
    return (original_value > threshold && original_value == dilated_value) ? 1u: 0u;
}

__device__ ushort4 is_peak(const float4 original_values, const float4 dilated_values)
{
    return ushort4 { (original_values.x == dilated_values.x) ? 1u: 0u,
                     (original_values.y == dilated_values.y) ? 1u: 0u,
                     (original_values.z == dilated_values.z) ? 1u: 0u,
                     (original_values.w == dilated_values.w) ? 1u: 0u };
}

__device__ float3 load_row(const Npp32f * input, int idx) {
    float3 val;
    val.x = __ldg(input+idx);
    val.y = __ldg(input+idx+1);
    val.z = __ldg(input+idx+2);
    return val;
}

__device__ float2 load_row2(const Npp32f * input, int idx) {
    float2 val = { __ldg(input+idx), __ldg(input+idx+1)};
    return val;
}


/**
 * loads a block of data. At the borders it loads values of 0.0f
 */
__device__ float3x3 load_block(const Npp32f * input, int idxX, int idxY, int width) {
    auto idx = idxY*width + idxX;
    float3x3 data;
    auto yrow = load_row(input, idx-1);
    auto xrow = load_row(input, idx-width-1);
    auto zrow = load_row(input, idx+width-1);
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

__device__ float3x3 load_block_top(const Npp32f * input, int idxX, int idxY, int width) {
     auto idx = idxY*width + idxX;
     float3x3 data;
     auto yrow = load_row(input, idx-1);
     auto zrow = load_row(input, idx+width-1);
     data.y1 = yrow.x;
     data.y2 = yrow.y;
     data.y3 = yrow.z;
     data.z1 = zrow.x;
     data.z2 = zrow.y;
     data.z3 = zrow.z;
     return data;
}

__device__ float3x3 load_block_bottom(const Npp32f * input, int idxX, int idxY, int width) {
     auto idx = idxY*width + idxX;
     float3x3 data;
     auto yrow = load_row(input, idx-1);
     auto xrow = load_row(input, idx-width-1);
     data.y1 = yrow.x;
     data.y2 = yrow.y;
     data.y3 = yrow.z;
     data.x1 = xrow.x;
     data.x2 = xrow.y;
     data.x3 = xrow.z;
     return data;
}

__device__ float3x3 load_block_left(const Npp32f * input, int idxX, int idxY, int width) {
     auto idx = idxY*width + idxX;
     float3x3 data;
     auto xrow = load_row2(input, idx-width);
     auto yrow = load_row2(input, idx);
     auto zrow = load_row2(input, idx+width);
     data.x2 = xrow.x;
     data.x3 = xrow.y;
     data.y2 = yrow.x;
     data.y3 = yrow.y;
     data.z2 = zrow.x;
     data.z3 = zrow.y;
     return data;
}

__device__ float3x3 load_block_right(const Npp32f * input, int idxX, int idxY, int width) {
     auto idx = idxY*width + idxX;
     float3x3 data;
     auto xrow = load_row2(input, idx-width-1);
     auto yrow = load_row2(input, idx-1);
     auto zrow = load_row2(input, idx+width-1);
     data.x1 = xrow.x;
     data.x2 = xrow.y;
     data.y1 = yrow.x;
     data.y2 = yrow.y;
     data.z1 = zrow.x;
     data.z2 = zrow.y;
     return data;
}

/**
 * Loads a 2x2 block starting at input and one row below, where the row width is specified by @p width
 */
__device__ float4 load_block_2x2(const Npp32f * input, int width) {
     float2 first = {__ldg(input), __ldg(input+1)};
     float2 second = {__ldg(input+width), __ldg(input+width+1)};
     return float4{first.x, first.y, second.x, second.y};
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

__global__ void find_peak(const Npp32f *d_input, const Npp32f * d_dilated, Npp16u * d_output)
{
    auto idx = blockDim.x * blockIdx.x + threadIdx.x;
    auto orig_values = reinterpret_cast<const float4*>(d_input)[idx];
    auto dilated_values = reinterpret_cast<const float4*>(d_dilated)[idx];
    (reinterpret_cast<ushort4*>(d_output))[idx] = is_peak(orig_values, dilated_values); 
}

/**
 * This version runs a thread per pixel
 * and hence needs to handle the special cases of the edges and the corners
 * where less data is loaded
 */
__global__ void dilate_peak_find(const Npp32f *d_input, Npp16u* d_output, const int width, const int height)
{
    auto idxX = blockDim.x * blockIdx.x + threadIdx.x;
    auto idxY = blockDim.y * blockIdx.y + threadIdx.y;
    if (idxX >= width) return;
    if (idxY >= height) return;

    unsigned short peak = 0u;
    //handle boundary conditions - top edge
    if (idxY == 0) {
        //Special case for height of 1 or width of 1
        if (width == 1 && height == 1) {
            peak = is_peak(d_input[0], d_input[0]);
        }
        //Top left corner case
        else if (idxX == 0) {
	    auto block = load_block_2x2(d_input, width);
	    auto dilated_value = dilate4(block);
	    peak = is_peak(block.x, dilated_value);
        } 
        //Top right corner case
        else if (idxX == (width-1)) {
	    auto block = load_block_2x2(d_input+width-2, width);
	    auto dilated_value = dilate4(block);
	    peak = is_peak(block.y, dilated_value);
	} else {
	    auto block = load_block_top(d_input, idxX, idxY, width);
            auto dilated_value = dilate3x3_top(block);
            peak = is_peak(block.y2, dilated_value);
        }
    //bottom edge
    } else if (idxY == height-1) {
        //Bottom left corner
        if (idxX == 0) {
	    auto block = load_block_2x2(d_input+width*(height-2), width);
            auto dilated_value = dilate4(block);
            peak = is_peak(block.z, dilated_value);
        }
        //Bottom right corner
        else if (idxX == (width-1)) {
	    auto block = load_block_2x2(d_input+width*(height-1)-2, width);
	    auto dilated_value = dilate4(block);
	    peak = is_peak(block.w, dilated_value);
        } else {
            auto block = load_block_bottom(d_input, idxX, idxY, width);        
            auto dilated_value = dilate3x3_bottom(block);
            peak = is_peak(block.y2, dilated_value);
        }
    //Left edge
    } else if (idxX == 0) {
        auto block = load_block_left(d_input, idxX, idxY, width);        
        auto dilated_value = dilate3x3_left(block);
        peak = is_peak(block.y2, dilated_value);
    
    //right edge
    } else if (idxX == (width-1)) {
        auto block = load_block_right(d_input, idxX, idxY, width);        
        auto dilated_value = dilate3x3_right(block);
        peak = is_peak(block.y2, dilated_value);

    } else {
        auto block = load_block(d_input, idxX, idxY, width);
        auto dilated_value = dilate3x3(block);
        peak = is_peak(block.y2, dilated_value);
    }
    d_output[idxY*width+idxX] = peak;
}

#define DILATE_BLOCK_SIZE 1024 
/**
 * This version uses a warp per block per row i.e.
 * 
 *  1 1 1 1 1 1 1 ... 1 2 2 2 2 2 ... 2
 *  3 3 3 3 3 3 3 ... 3 4 4 4 4 4 ... 4
 *  5 5 5 5 5 5 5 ... 5 6 6 6 6 6 ... 6
 *
 * where the number represents the block index.
 * The warp collaboratively loads the data for each row
 * and then dilates each pixel and then writes out the
 * results in a co-operative manner.
 */
template<int dilate_block_size=DILATE_BLOCK_SIZE>
__global__ void dilate_peak_find_v2(const Npp32f * d_input, Npp16u* d_output, const int width, const float threshold)
{
    //TODO: try swapping block x and y - however order of block execution is undefined so this may not improve cache behaviour - hence v3 impl to solve this
    auto idxX = threadIdx.x; // Index into shared memory cache
    auto gidxX = blockDim.x * blockIdx.x + idxX;

    if (gidxX >= width) return;

    auto row = blockIdx.y;
    bool is_last_row = row == (gridDim.y-1);
    bool is_first_row = (row == 0);
    bool is_leftmost_element = (idxX == 0);
    bool is_rightmost_element = ((idxX == (blockDim.x-1)) || idxX == (width-1));
    bool is_first_block = (blockIdx.x == 0);
    bool is_last_block = (blockIdx.x == (gridDim.x-1));

    __shared__ float data[3][dilate_block_size];
    //Always Load the middle row
    data[1][idxX] = d_input[row*width+gidxX];
    if (!is_first_row) {
        data[0][idxX] = d_input[(row-1)*width+gidxX];
    }
    if (!is_last_row) {
        data[2][idxX] = d_input[(row+1)*width+gidxX];
    }
    
    auto dilated_value = data[1][idxX];
    if (!is_first_row) {
        dilated_value = fmaxf(dilated_value, data[0][idxX]);
    }
    if (!is_last_row) {
        dilated_value = fmaxf(dilated_value, data[2][idxX]);
    }

    //Share the dilated values via shared memory
    data[0][idxX] = dilated_value;

    __syncthreads();

    //Left hand boundary condition
    float cmp_val_l = 0.0f;
    float cmp_val_r = 0.0f;
    if (is_leftmost_element) {
        if (not is_first_block) {
	    float temp1 = 0.0f; 
	    float temp2 = 0.0f; 
	    if (not is_first_row) {
                temp1 = __ldg(d_input+(row-1)*width+gidxX-1);
	    }
	    if (not is_last_row) {
                temp2 = __ldg(d_input+(row+1)*width+gidxX-1);
	    }
       	    cmp_val_l = fmaxf(temp1, temp2); 
	    cmp_val_l = fmaxf(cmp_val_l, __ldg(d_input+row*width+gidxX-1));
        }
    } else {
        cmp_val_l = data[0][idxX-1];
    }
    //Right hand boundary condition
    if (is_rightmost_element) {
        if (not is_last_block) {
	    float temp1 = 0.0f; 
	    float temp2 = 0.0f; 
	    if (not is_first_row) {
                temp1 = __ldg(d_input+(row-1)*width+gidxX+1);
	    }
	    if (not is_last_row) {
                temp2 = __ldg(d_input+(row+1)*width+gidxX+1);
	    }
       	    cmp_val_r = fmaxf(temp1, temp2); 
	    cmp_val_r = fmaxf(cmp_val_r, __ldg(d_input+row*width+gidxX-1));
        }
    } else {
        cmp_val_r = data[0][idxX+1];
    }
    dilated_value = fmaxf(dilated_value, cmp_val_l);
    dilated_value = fmaxf(dilated_value, cmp_val_r);

    d_output[row*width+gidxX] = is_peak(data[1][idxX], dilated_value, threshold);
}

/**
 * This version uses a block of warps laid out as follows: 
 * 
 *   1  2  3 .... 32
 *  33 34 35 .... 64 
 *  ..       .... .. 
 *  ..       .... 1024 
 *
 * where the number represents the thread index.
 * The warp collaboratively loads the data for each row
 * and then dilates each pixel and then writes out the
 * results in a co-operative manner.
 * This is better than v2 because it facilitates data
 * reuse (each warp only loads a single row - except the first and last)
 * and the left and right boundary conditions are no longer divergent 
 */
__global__ void dilate_peak_find_v3(const Npp32f * d_input, Npp16u* d_output, const int linestep, const int width, const int height, const float threshold)
{
    const auto idxX = threadIdx.x; // Index into shared memory cache
    const auto gidxX = blockIdx.x * blockDim.x + idxX;
    const auto row_s = threadIdx.y + 1;
    const auto row = (blockIdx.y * blockDim.y + threadIdx.y);
    const bool is_last_row_in_block = (row_s == blockDim.y);
    const bool is_last_row = row == (height-1); 
    const bool is_last_active_row_in_block = is_last_row ? true: is_last_row_in_block;
    const bool is_first_row = (row == 0);
    const bool is_first_row_in_block = (row_s == 1);
    const bool is_leftmost_element = (gidxX == 0);
    const bool is_leftmost_element_in_block = (idxX == 0);
    const bool is_rightmost_element = gidxX == (width-1);
    const bool is_rightmost_element_in_block = ((idxX == blockDim.x-1) || is_rightmost_element);
    const int linewidth = width+linestep;
    
    __shared__ float data[34][32];
    float dilated_value, my_val = 0.0f;
    float val1 = 0.0f; 
    float val2 = 0.0f; 
 
    const bool thread_active = !(gidxX >= width || row >= height); 
    const auto last_row_s = (blockIdx.y == (gridDim.y-1) ? (height-1) & (blockDim.y-1) : blockDim.y + 1);
    const auto idx_s = (is_first_row_in_block) ? 0 : last_row_s;

    //Note we can't return here as we use syncthreads
    if (thread_active) {

        //Load block data collaboratively
        auto * ptr = d_input + row*linewidth+gidxX;
        my_val = __ldg(ptr);

	//Boundary condition neightbours on my row
        if ((is_leftmost_element_in_block and not is_leftmost_element)
           || (is_rightmost_element_in_block and not is_rightmost_element))
            val2 = __ldg( is_leftmost_element_in_block ? ptr-1 : ptr+1);

	//Top/bottom edge boundary condition
        if ((is_first_row_in_block and not is_first_row)
           || (is_last_active_row_in_block and not is_last_row)) {
            auto ptr_bc = is_first_row_in_block ? ptr - linewidth : ptr + linewidth;
	    data[idx_s][idxX] = __ldg(ptr_bc);
            if ((is_leftmost_element_in_block and not is_leftmost_element)
               || (is_rightmost_element_in_block and not is_rightmost_element)) {
                val1 = __ldg(is_leftmost_element_in_block ? ptr_bc - 1: ptr_bc + 1);
            }
        }

        //Compare against left and right values
        auto tmp = fmaxf(__shfl_up(my_val, 1), __shfl_down(my_val, 1));

	//Compare against boundary condition value
        auto tmp2 = fmaxf(my_val, val2);
        dilated_value = fmaxf(tmp, tmp2);

        //Share the dilated values via shared memory
        data[row_s][idxX] = dilated_value;

        // Dilate the first and last rows
	if ((is_first_row_in_block and not is_first_row)
           || (is_last_active_row_in_block and not is_last_row)) {
            float lval;
            auto dilated_value_bc = lval = data[idx_s][idxX];
            dilated_value_bc = fmaxf(dilated_value_bc, __shfl_down(lval, 1));
            dilated_value_bc = fmaxf(dilated_value_bc, __shfl_up(lval, 1));
            dilated_value_bc = fmaxf(dilated_value_bc, val1);
            data[idx_s][idxX] = dilated_value_bc;
        }
    }
    //Sync here so we can see the data from other warps i.e. the other rows in this block
    __syncthreads();
 
    if (thread_active) {
        if (not is_first_row) {
            dilated_value = fmaxf(dilated_value, data[row_s-1][idxX]);
        }
        if (not is_last_row) {
            dilated_value = fmaxf(dilated_value, data[row_s+1][idxX]);
        }

        d_output[row*linewidth+gidxX] = is_peak(my_val, dilated_value, threshold);
    }
}

/**
 * V4 runs in the same manner as v3 but avoids the top bottom
 * boundary condition by running more warps
 * Each block now processes 30 rows instead of 32 and the 0 and 31st warp
 * 'handle' the boundary condition 
 */
__global__ void dilate_peak_find_v4(const Npp32f * d_input, Npp16u* d_output, const int linestep, const int width, const int height, const float threshold)
{
    const auto idxX = threadIdx.x; // Index into shared memory cache
    const auto gidxX = blockIdx.x * blockDim.x + idxX;
    const auto row_s = threadIdx.y;
    const auto row = ((blockIdx.y) * (blockDim.y-2)) + threadIdx.y;
    const bool is_last_row = row == (height-1); 
    const bool is_first_row = (row == 0);
    const bool is_leftmost_element = (gidxX == 0);
    const bool is_leftmost_element_in_block = (idxX == 0);
    const bool is_rightmost_element = gidxX == (width-1);
    const bool is_rightmost_element_in_block = ((idxX == blockDim.x-1) || is_rightmost_element);
    const int linewidth = width+linestep;
    
    __shared__ float data[32][32];
    float dilated_value, my_val = 0.0f;
    float val1 = 0.0f; 
 
    const bool thread_active = !(gidxX >= width || row >= height); 
    const bool thread_output = thread_active && ((threadIdx.y > 0 && threadIdx.y < 31) || is_first_row || is_last_row);

    //Note we can't return here as we use syncthreads
    if (thread_active) {

        //Load block data collaboratively
        auto * ptr = d_input + row*linewidth + gidxX;
        my_val = __ldg(ptr);

	//Boundary condition neightbours on my row
        if ((is_leftmost_element_in_block and not is_leftmost_element)
           || (is_rightmost_element_in_block and not is_rightmost_element))
            val1 = __ldg( is_leftmost_element_in_block ? ptr-1 : ptr+1);

        //Compare against left and right values
        auto tmp = fmaxf(__shfl_up(my_val, 1), __shfl_down(my_val, 1));

	//Compare against boundary condition value
        auto tmp2 = fmaxf(my_val, val1);
        dilated_value = fmaxf(tmp, tmp2);

        //Share the dilated values via shared memory
        data[row_s][idxX] = dilated_value;
    }
    //Sync here so we can see the data from other warps i.e. the other rows in this block
    __syncthreads();
 
    if (thread_output) {
        if ( not is_first_row ) {
            dilated_value = fmaxf(dilated_value, data[row_s-1][idxX]);
        }
        if ( not is_last_row ) {
            dilated_value = fmaxf(dilated_value, data[row_s+1][idxX]);
        }

        d_output[row*linewidth+gidxX] = is_peak(my_val, dilated_value, threshold);
    }
}


class PeakFinderContext
{
public:
    PeakFinderContext(int width, int height)
	//: dilated(width, height)
    {
	cudaStreamCreate(&stream);
    }

    PeakFinderContext(const PeakFinderContext &) = delete;
    PeakFinderContext & operator = (const PeakFinderContext &) = delete;

    ~PeakFinderContext() {
	cudaStreamDestroy(stream);
    } 

    void v1(const NppImage & input, unsigned short * output) {
	runFusedKernel(input, output);
	cudaStreamSynchronize(stream);
    }

    void v2(const NppImage & input, unsigned short * output, const float threshold) {
#if 0
	int blockSize;   // The launch configurator returned block size 
	int minGridSize; // The minimum grid size needed to achieve the 
	// maximum occupancy for a full device launch 

	// Attach the used memory to this stream
	cudaStreamAttachMemAsync(this->stream, input.data);
	cudaStreamAttachMemAsync(this->stream, output);

        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, 
                                           find_peak, 0, 0); 
	// Round up according to array size 
	int gridSize = ((input.size().width * input.size().height) + blockSize - 1) / blockSize; 
	gridSize /= 4; //we operate on four elements at a time
        std::cout << "Gridsize " << gridSize << std::endl;	
	runKernels(input, output, gridSize, blockSize, this->stream);
#endif
	runFusedKernelv2(input, output, threshold);
	cudaStreamSynchronize(stream);
    }

    void v3(const NppImage & input, unsigned short * output, const float threshold) {
        dim3 blockDim = {32, 32, 1};
        dim3 gridSize = {1u + ((input.size().width-1)/blockDim.x), 1u + ((input.size().height-1)/blockDim.y), 1};
	//std::cout << "Grid dims: " << gridSize.x << " " << gridSize.y << " " << gridSize.z << std::endl;	
	dilate_peak_find_v3<<<gridSize, blockDim, 0, stream>>>(input.data, output, input.linestep, input.size().width, input.size().height, threshold);
 	cudaStreamSynchronize(stream);
    }

    void v4(const NppImage & input, unsigned short * output, const float threshold) {
        dim3 blockDim = {32, 32, 1};
        dim3 gridSize = {1u + ((input.size().width-1)/blockDim.x), 1u + ((input.size().height-1)/(blockDim.y-2)), 1};
	//std::cout << "Grid dims: " << gridSize.x << " " << gridSize.y << " " << gridSize.z << std::endl;
	dilate_peak_find_v4<<<gridSize, blockDim, 0, stream>>>(input.data, output, input.linestep, input.size().width, input.size().height, threshold);
 	cudaStreamSynchronize(stream);
    }


private:

    void runFusedKernelv2(const NppImage & input, unsigned short * output, const float threshold) {
        dim3 blockDim = {DILATE_BLOCK_SIZE, 1, 1};
        dim3 gridSize = {1 + ((input.size().width-1)/blockDim.x), 1 + ((input.size().height-1)/blockDim.y), 1};
	//std::cout << "Grid dims: " << gridSize.x << " " << gridSize.y << " " << gridSize.z << std::endl;	
	dilate_peak_find_v2<><<<gridSize, blockDim, 0, stream>>>(input.data, output, input.size().width, threshold);
    }

    void runFusedKernel(const NppImage & input, unsigned short * output) {
	int blockSize;   // The launch configurator returned block size 
	int minGridSize; // The minimum grid size needed to achieve the 
	// maximum occupancy for a full device launch 

	// Attach the used memory to this stream
	cudaStreamAttachMemAsync(stream, input.data);
	cudaStreamAttachMemAsync(stream, output);

        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, 
                                           dilate_peak_find, 0, 0);

	dim3 blockDim = {32, 2, 1};
        dim3 gridSize = {1 + ((input.size().width-1)/blockDim.x), 1 + ((input.size().height-1)/blockDim.y), 1};
        dilate_peak_find<<<gridSize, blockDim, 0, stream>>>(input.data, output, input.size().width, input.size().height);
    }

    void runKernels(const NppImage & input, unsigned short * output, int gridSize, int blockSize, cudaStream_t stream) {
	//std::cout << "Gridsize: " << gridSize << " width " << input.size().width << " height " << input.size().height << " blockSize " << blockSize << std::endl;

        AllocatedNppImage dilated(input.size().width, input.size().height);
        NppiPoint srcOffset = {0, 0};
        //Do a 3x3 image dilation
	nppSetStream(stream);
	auto status = nppiDilate3x3Border_32f_C1R(input.data, input.linestep, input.size(), srcOffset, dilated.data, dilated.linestep, dilated.size(), NPP_BORDER_REPLICATE);
        if (status != NPP_SUCCESS) handle_npp_error(status);
	find_peak<<<gridSize, blockSize, 0, stream>>>(input.data, dilated.data, output);
    }

    cudaStream_t stream;
};

void peakfind_v3(float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output)
{
   	if (height == 0 || width == 0) return;	
	NppImage input(width, height, d_input, d_input_linestep);
	PeakFinderContext p(width, height);
        p.v3(input, d_output, 0.0f);
}

void peakfind(float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output)
{
   	if (height == 0 || width == 0) return;	
	NppImage input(width, height, d_input, d_input_linestep);
	PeakFinderContext p(width, height);
        p.v1(input, d_output);
}

void peakfind_v2(float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output)
{
   	if (height == 0 || width == 0) return;	
	NppImage input(width, height, d_input, d_input_linestep);
	PeakFinderContext p(width, height);
        p.v2(input, d_output, 0.0f);
}

