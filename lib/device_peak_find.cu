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

__global__ void find_peak(const Npp32f *d_input, const Npp32f * d_dilated, Npp16u * d_output)
{
    auto idx = blockDim.x*blockIdx.x + threadIdx.x;
/*    auto orig_values = __ldg(reinterpret_cast<const float4*>(d_input) + idx);
    auto dilated_values = __ldg(reinterpret_cast<const float4*>(d_dilated) + idx);*/
    auto orig_values = reinterpret_cast<const float4*>(d_input)[idx];
    auto dilated_values = reinterpret_cast<const float4*>(d_dilated)[idx];
    unsigned short is_peak  = (orig_values.x == dilated_values.x) ? 1u : 0u;
    unsigned short is_peak2 = (orig_values.y == dilated_values.y) ? 1u : 0u;
    unsigned short is_peak3  = (orig_values.z == dilated_values.z) ? 1u : 0u;
    unsigned short is_peak4 = (orig_values.w == dilated_values.w) ? 1u : 0u;
    (reinterpret_cast<ushort4*>(d_output))[idx] = ushort4{is_peak, is_peak2, is_peak3, is_peak4};
}

class PeakFinderContext
{
public:
    PeakFinderContext(int width, int height)
	: dilated(width, height)
    {
	cudaStreamCreate(&stream);
    }

    PeakFinderContext(const PeakFinderContext &) = delete;
    PeakFinderContext & operator = (const PeakFinderContext &) = delete;

    ~PeakFinderContext() {
	cudaStreamDestroy(stream);
    } 

    void operator()(const NppImage & input, unsigned short * output) {
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

    void runKernels(const NppImage & input, unsigned short * output, int gridSize, int blockSize, cudaStream_t stream) {
	//std::cout << "Gridsize: " << gridSize << " width " << input.size().width << " height " << input.size().height << " blockSize " << blockSize << std::endl;
        NppiPoint srcOffset = {0, 0};
        //Do a 3x3 image dilation
	nppSetStream(stream);
	auto status = nppiDilate3x3Border_32f_C1R(input.data, input.linestep, input.size(), srcOffset, dilated.data, dilated.linestep, dilated.size(), NPP_BORDER_REPLICATE);
        if (status != NPP_SUCCESS) handle_npp_error(status);
	find_peak<<<gridSize, blockSize, 0, stream>>>(input.data, dilated.data, output);
    }

    AllocatedNppImage dilated;
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

