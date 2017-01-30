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

__device__ Npp16u is_peak(const float original_value, const float dilated_value)
{
    return (original_value == dilated_value) ? 1u: 0u;
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
    float2 val = __ldg((float2*)(input+idx));
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
     float2 first = __ldg((float2 *)input);
     float2 second = __ldg((float2 *)(input+width));
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

__global__ void dilate_peak_find(const Npp32f *d_input, Npp16u* d_output)
{
    auto idxX = blockDim.x * blockIdx.x + threadIdx.x;
    auto idxY = blockDim.y * blockIdx.y + threadIdx.y;
    const int width = blockDim.x * gridDim.x;
    //handle boundary conditions - top edge
    if (blockIdx.y == 0 && threadIdx.y == 0) {
        //Top left corner case
        if (blockIdx.x == 0 && threadIdx.x == 0) {
	    auto block = load_block_2x2(d_input, width);
	    auto dilated_value = dilate4(block);
	    auto peak = is_peak(block.x, dilated_value);
	    d_output[idxY*width+idxX] = peak;
        } 
        //Top right corner case
        else if (blockIdx.x == (gridDim.x-1) && threadIdx.x == (blockDim.x-1)) {
	    auto block = load_block_2x2(d_input+width-2, width);
	    auto dilated_value = dilate4(block);
	    auto peak = is_peak(block.y, dilated_value);
	    d_output[idxY*width+idxX] = peak;
	} else {
	    auto block = load_block_top(d_input, idxX, idxY, width);
            auto dilated_value = dilate3x3_top(block);
            auto peak = is_peak(block.y2, dilated_value);
            d_output[idxY*width+idxX] = peak;
        }
    //bottom edge
    } else if (blockIdx.y == (gridDim.y-1) && threadIdx.y == blockDim.y-1) {
        //Bottom left corner
        if (blockIdx.x == 0 && threadIdx.x == 0) {
	    auto block = load_block_2x2(d_input+(width*((blockDim.y-1)*(gridDim.y-1))), width);
            auto dilated_value = dilate4(block);
            auto peak = is_peak(block.z, dilated_value);
            d_output[idxY*width+idxX] = peak;
        }
        //Bottom right corner
        else if (blockIdx.x == (gridDim.x-1) && threadIdx.x == blockDim.x-1) {
	    auto block = load_block_2x2(d_input+(width*(1+(blockDim.y-1)*(gridDim.y-1)))-2, width);
	    auto dilated_value = dilate4(block);
	    auto peak = is_peak(block.w, dilated_value);
	    d_output[idxY*width+idxX] = peak;
        } else {
            auto block = load_block_bottom(d_input, idxX, idxY, width);        
            auto dilated_value = dilate3x3_bottom(block);
            auto peak = is_peak(block.y2, dilated_value);
            d_output[idxY*width+idxX] = peak;
        }
    //Left edge
    } else if (blockIdx.x == 0 && threadIdx.x == 0) {
        auto block = load_block_left(d_input, idxX, idxY, width);        
        auto dilated_value = dilate3x3_left(block);
        auto peak = is_peak(block.y2, dilated_value);
        d_output[idxY*width+idxX] = peak;

    
    //right edge
    } else if (blockIdx.x == (gridDim.x-1) && threadIdx.x == blockDim.x-1) {
        auto block = load_block_right(d_input, idxX, idxY, width);        
        auto dilated_value = dilate3x3_right(block);
        auto peak = is_peak(block.y2, dilated_value);
        d_output[idxY*width+idxX] = peak;

     } else {
        auto block = load_block(d_input, idxX, idxY, width);
        auto dilated_value = dilate3x3(block);
        auto peak = is_peak(block.y2, dilated_value);
        d_output[idxY*width+idxX] = peak;
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

    void operator()(const NppImage & input, unsigned short * output) {
	runFusedKernel(input, output);
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
	runKernels(input, output, gridSize, blockSize, this->stream);
#endif
	cudaStreamSynchronize(stream);
    }

    void operator()(const NppImage & input, unsigned short * output, cudaStream_t stream) {
	int blockSize;   // The launch configurator returned block size 
	int minGridSize; // The minimum grid size needed to achieve the 
	// maximum occupancy for a full device launch 

	// Attach the used memory to this stream
	cudaStreamAttachMemAsync(stream, input.data);
	cudaStreamAttachMemAsync(stream, output);

        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, 
                                           find_peak, 0, 0); 
	// Round up according to array size 
	int gridSize = ((input.size().width * input.size().height) + blockSize - 1) / blockSize; 
	gridSize /= 4; //we operate on four elements at a time
	runKernels(input, output, gridSize, blockSize, stream);
    }

private:

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
	dim3 gridSize = {input.size().width/blockDim.x, std::max(1u, input.size().height/blockDim.y), 1};
        dilate_peak_find<<<gridSize, blockDim, 0, stream>>>(input.data, output);
    }

    void runKernels(const NppImage & input, unsigned short * output, int gridSize, int blockSize, cudaStream_t stream) {
	//std::cout << "Gridsize: " << gridSize << " width " << input.size().width << " height " << input.size().height << " blockSize " << blockSize << std::endl;
/*
        NppiPoint srcOffset = {0, 0};
        //Do a 3x3 image dilation
	nppSetStream(stream);
	auto status = nppiDilate3x3Border_32f_C1R(input.data, input.linestep, input.size(), srcOffset, dilated.data, dilated.linestep, dilated.size(), NPP_BORDER_REPLICATE);
        if (status != NPP_SUCCESS) handle_npp_error(status);
	find_peak<<<gridSize, blockSize, 0, stream>>>(input.data, dilated.data, output);
*/
    }

    //AllocatedNppImage dilated;
    cudaStream_t stream;
};

class PeakFinderManager
{
public:
	PeakFinderManager(int count, int width, int height)
	{
            for(int i=0; i < count; ++i) {
		peak_finders.emplace_back(new PeakFinderContext(width, height));
                free_peak_finders.push_back(peak_finders.back().get());
	    }
	}

	~PeakFinderManager() {
	    //Everyone should have released their PeakFinderContext's before
            // this is called
            assert(peak_finders.size() == free_peak_finders.size());
        }

	PeakFinderContext * getFinderContext() {
	    auto context = free_peak_finders.back();
            free_peak_finders.pop_back();
            return context;
        }

        void releaseFinderContext(PeakFinderContext * context) {
	    free_peak_finders.push_back(context);
        }
private:
	std::vector<std::unique_ptr<PeakFinderContext>> peak_finders;
	std::vector<PeakFinderContext *> free_peak_finders;
};


class PeakFinder
{
public:
	PeakFinder(PeakFinderManager & m) 
	: manager(m), context(m.getFinderContext())
	{
	}

	~PeakFinder()
	{
	    manager.releaseFinderContext(context);
	}

	PeakFinder(const PeakFinder &) = delete;
	PeakFinder & operator= (const PeakFinder &) = delete;

	void operator()(const NppImage & image, unsigned short * output)
	{
	    (*context)(image, output);
	}
private:
	PeakFinderManager & manager;
	PeakFinderContext * context;
};

void peakfind(float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output)
{
	/*PeakFinderManager peakFinderManager(1, 1024, 1024);
	PeakFinder p(peakFinderManager);*/
	NppImage input(width, height, d_input, d_input_linestep);
	PeakFinderContext p(width, height);
        p(input, d_output);
}

//This version isn't quite right as it doesn't keep the PeakFinderContext alive and so the dilated image will get deallocated too early
/*void peakfindAsync(float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output, cudaStream_t stream)
{
	NppImage input(width, height, d_input, d_input_linestep);
	PeakFinderContext p(width, height);
        p(input, d_output, stream);
}*/

