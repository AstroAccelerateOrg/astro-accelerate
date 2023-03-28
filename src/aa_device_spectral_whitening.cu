#include <stdio.h>
#include <cufft.h>
#include <vector>
#include "aa_params.hpp"
#include "aa_device_spectrum_whitening_kernel.hpp"
#include "aa_host_utilities.hpp"

#include "aa_log.hpp"

// define for debug info
//#define AA_SPECTRAL_WHITENING_DEBUG

//{{{ Dopler Stretch 

namespace astroaccelerate {

	void create_dered_segment_sizes_prefix_sum(std::vector<int> *segment_sizes, int min_segment_length, int max_segment_length, size_t nSamples){
		segment_sizes->push_back(1);
		
		int seg_length = min_segment_length;
		int seg_sum = 1 + seg_length;
		int new_seg_length = min_segment_length*log(seg_sum);
		if(new_seg_length > max_segment_length) new_seg_length = max_segment_length;
		segment_sizes->push_back(seg_sum);
		
		while( (size_t) (seg_sum + seg_length) < nSamples) {
			seg_length = new_seg_length;
			seg_sum = seg_sum + seg_length;
			new_seg_length = min_segment_length*log(seg_sum);
			if(new_seg_length > max_segment_length) new_seg_length = max_segment_length;
			segment_sizes->push_back(seg_sum);
		}
		segment_sizes->push_back(nSamples);
	}
	
	void calculate_mean_for_segments(float *d_MSD_segments, float2 *d_input, size_t nSamples, int nDMs, int *d_segment_sizes, int nSegments, cudaStream_t &stream){
		size_t smem_size = 0;
		dim3 grid_size(nSegments, nDMs , 1);
		dim3 block_size(256, 1, 1);
		call_kernel_segmented_MSD(block_size, grid_size, smem_size, stream, d_MSD_segments, d_input, d_segment_sizes, nSamples, nSegments);
	}
	
	void calculate_median_for_segments(float *d_MSD_segments, float2 *d_input, size_t nSamples, int nDMs, int *d_segment_sizes, int nSegments, cudaStream_t &stream){
		size_t smem_size = 0;
		dim3 grid_size(nSegments - 1, nDMs , 1);
		dim3 block_size(128, 1, 1);
		call_kernel_segmented_median(block_size, grid_size, smem_size, stream, d_MSD_segments, d_input, d_segment_sizes, nSamples, nSegments);
	}

	void spectrum_whitening_SGP1(float *d_input, size_t nSamples, int nDMs, cudaStream_t &stream){
		dim3 grid_size((nSamples + 128 - 1)/128, nDMs , 1);
		dim3 block_size(128, 1, 1);
		call_kernel_spectrum_whitening_SGP1(
			block_size, 
			grid_size, 
			0, 
			stream, 
			d_input, 
			nSamples
		);
	}
	
	void spectrum_whitening_SGP2(float2 *d_input, size_t nSamples, int nDMs, bool enable_median, cudaStream_t &stream) {
		int max_segment_length = 200;
		int min_segment_length = 6;
		
		std::vector<int> segment_sizes; // must be a prefix sum!
		create_dered_segment_sizes_prefix_sum(&segment_sizes, min_segment_length, max_segment_length, nSamples);
		
		//------------ Allocate and copy segment sizes to the GPU
		int nSegments = segment_sizes.size();
		size_t ss_size = nSegments*sizeof(int);
		size_t ss_MSD_size = nSegments*sizeof(float)*nDMs;
		int *d_segment_sizes;
		if ( cudaSuccess != cudaMalloc((void **) &d_segment_sizes, ss_size )) {
			printf("spectral whitening Allocation error! d_segment_sizes\n");
		}
		float *d_MSD_segments;
		if ( cudaSuccess != cudaMalloc((void **) &d_MSD_segments, ss_MSD_size )) {
			printf("spectral whitening Allocation error! d_MSD_segments\n");
		}
		cudaError_t e = cudaMemcpy(d_segment_sizes, segment_sizes.data(), ss_size, cudaMemcpyHostToDevice);
		if(e != cudaSuccess) {
			LOG(log_level::error, "Could not cudaMemcpy in aa_device_spectral_whitening.cu (" + std::string(cudaGetErrorString(e)) + ")");
		}
		
		//------------ Calculate mean or median;
		if(enable_median){
			calculate_median_for_segments(d_MSD_segments, d_input, nSamples, nDMs, d_segment_sizes, nSegments, stream);
		}
		else {
			calculate_mean_for_segments(d_MSD_segments, d_input, nSamples, nDMs, d_segment_sizes, nSegments, stream);
		}
		
		#ifdef AA_SPECTRAL_WHITENING_DEBUG
		float *h_segmented_MSD;
		h_segmented_MSD = new float[nDMs*nSegments];
		e = cudaMemcpy(h_segmented_MSD, d_MSD_segments, nDMs*nSegments*sizeof(float), cudaMemcpyDeviceToHost);
		printf("Calculated means:\n");
		char str[300];
		for(int d=0; d<nDMs; d++){
			if(enable_median){
				sprintf(str, "PSR_fft_GPU_medians_%d.dat", d);
			}
			else {
				sprintf(str, "PSR_fft_GPU_means_%d.dat", d);
			}
			Export_data_to_file(&h_segmented_MSD[d*nSegments], nSegments, 1, str);
		}
		delete[] h_segmented_MSD;
		#endif
		
		//------------ Call kernel for spectral whitening
		size_t smem_size = 0;
		dim3 grid_size(nSegments-1, nDMs , 1);
		dim3 block_size(256, 1, 1);
		call_kernel_spectrum_whitening_SGP2(block_size, grid_size, smem_size, stream, d_MSD_segments, d_input, d_segment_sizes, nSamples, nSegments);
		
		cudaFree(d_segment_sizes);
		cudaFree(d_MSD_segments);
	}

} //namespace astroaccelerate
