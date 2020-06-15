#ifndef ASTRO_ACCELERATE_AA_JERK_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_JERK_STRATEGY_HPP

#include "aa_strategy.hpp"
#include "aa_jerk_plan.h"
#include <iostream>
#include <vector>

// NOTES:
// conv_size should really be determined by the class based on size of the filter and performance of the GPU
namespace astroaccelerate {
	
/** \class aa_jerk_strategy aa_jerk_strategy.hpp "include/aa_jerk_strategy.hpp"
 * \brief Class to configure an jerk search strategy.
 * \details An jerk strategy is required for running any pipeline that will run the fdas/jerk search component aa_pipeline::component::jerk.
 * \details This class calculates memory requirements and configures the jerk search component.
 * \author AstroAccelerate.
 * \date 2020-06-15.
 */ 
class aa_jerk_strategy : public aa_strategy {
private:
	//-----------------> JERK plan
	// input parameters
	size_t nSamples_time_dom;
	size_t nSamples_freq_dom;
	size_t nSamples_output_plane;
	size_t nDMs;
	
	// filter parameters
	float z_max_search_limit;
	float z_search_step;
	float w_max_search_limit;
	float w_search_step;
	int interbinned_samples;
	int high_precision;
	
	// output parameters
	bool MSD_outlier_rejection;
	int candidate_selection;
	
	//-----------------> JERK strategy
	
	int conv_size; // size of the segment in overlap and save
	int nFilters_z_half;
	int nFilters_z;
	int nFilters_w_half;
	int nFilters_w;
	int nFilters_total;
	
	int filter_length;
	int filter_halfwidth;
	int useful_part_size;
	
	int nSegments;
	size_t output_size_one_DM;
	size_t output_size_z_plane;
	size_t output_size_total;
	
	int nZPlanes;
	int nZPlanes_per_chunk;
	int nZPlanes_chunks;
	int nZPlanes_remainder;
	std::vector<int> ZW_chunks;
	size_t free_memory;
	size_t required_memory;
	size_t total_memory;
	size_t reserved_memory_for_candidate_selection;
	
	size_t filter_padded_size;
	size_t filter_padded_size_bytes;
	
	bool c_ready;
	
	void calculate_memory_split(size_t available_free_memory){
		//---------> Memory testing and splitting
		nZPlanes = nFilters_w;
		size_t z_plane_size_bytes = output_size_z_plane*sizeof(float);
		available_free_memory = available_free_memory - reserved_memory_for_candidate_selection - filter_padded_size_bytes;
		nZPlanes_per_chunk = (int) (available_free_memory/z_plane_size_bytes);
		nZPlanes_chunks = (int) (nZPlanes/nZPlanes_per_chunk);
		nZPlanes_remainder = nZPlanes - nZPlanes_chunks*nZPlanes_per_chunk;
		
		for(int f=0; f<nZPlanes_chunks; f++){
			ZW_chunks.push_back(nZPlanes_per_chunk);
		}
		if(nZPlanes_remainder>0) {
			ZW_chunks.push_back(nZPlanes_remainder);
		}
		
		required_memory = (nZPlanes_chunks + nZPlanes_remainder)*z_plane_size_bytes + reserved_memory_for_candidate_selection;
	}
	
public:
	aa_jerk_strategy(aa_jerk_plan plan, size_t available_free_memory){
		nSamples_time_dom   = plan.nTimesamples;
		nSamples_freq_dom   = (plan.nTimesamples>>1) + 1; //because R2C FFT
		nDMs                = plan.nDMs;
		
		// Calculate number of filters
		// number of filters must also account for negative accelerations and w=z=0;
		nFilters_z_half   = plan.z_max_search_limit/plan.z_search_step;
		nFilters_z        = nFilters_z_half + nFilters_z_half + 1; 
		nFilters_w_half   = plan.w_max_search_limit/plan.w_search_step;
		nFilters_w        = nFilters_w_half + nFilters_w_half + 1;
		nFilters_total    = nFilters_z*nFilters_w;
		
		// recompute maximum z and w values based on step
		z_max_search_limit  = nFilters_z_half*plan.z_search_step;
		z_search_step       = plan.z_search_step;
		w_max_search_limit  = nFilters_w_half*plan.w_search_step;
		w_search_step       = plan.w_search_step;
		interbinned_samples = plan.interbinned_samples;
		high_precision      = plan.high_precision;
		
		// Strategy
		conv_size = 2048;
		
		filter_length    = plan.filter_length*interbinned_samples;
		filter_halfwidth = plan.filter_halfwidth*interbinned_samples;
		useful_part_size = conv_size - 2*filter_halfwidth + 1;
		
		nSegments           = (nSamples_freq_dom + useful_part_size - 1)/useful_part_size;
		output_size_one_DM  = nSegments*useful_part_size;
		output_size_z_plane = nFilters_z*output_size_one_DM;
		output_size_total   = nFilters_total*output_size_one_DM;
		reserved_memory_for_candidate_selection = 2*output_size_z_plane*sizeof(float);
		
		filter_padded_size       = nFilters_total*conv_size;
		filter_padded_size_bytes = nFilters_total*conv_size*sizeof(float2); //*8 for complex float
		
		calculate_memory_split(available_free_memory);
	}
	
	~aa_jerk_strategy(){
		ZW_chunks.clear();
	}
	
    bool ready() const {
      return c_ready;
    }
	
    bool setup() {
      return ready();
    }
	
	void PrintStrategy(){
		printf("-------------------------------------------\n");
		printf("Input parameters:\n");
		printf("    Number of time samples before FFT: %zu\n", nSamples_time_dom);
		printf("    Number of time samples after FFT:  %zu\n", nSamples_freq_dom);
		printf("    Number of time samples in output:  %zu\n", nSamples_freq_dom);
		printf("    Number of DM trials:               %zu\n", nDMs);
		printf("Filter parameters:\n");
		printf("    Filter's length    %d;\n", filter_length);
		printf("    Filter's halfwidth %d;\n", filter_halfwidth);
		printf("    z max:             %f;\n", z_max_search_limit);
		printf("    z step size:       %f;\n", z_search_step);
		printf("    w max:             %f;\n", w_max_search_limit);
		printf("    w step size:       %f;\n", w_search_step);
		printf("\n");
		printf("Interbinning: ");
		if(interbinned_samples==2) printf("Yes.\n"); else printf("No.\n");
		printf("High precision filters: ");
		if(high_precision==1) printf("Yes.\n"); else printf("No.\n");
		printf("-------------------------------------------\n");
		printf("\n");
		printf("-------------------------------------------\n");
		printf("Convolution size: %d\n", conv_size);
		printf("Half filters widths z=%d; w=%d\n", nFilters_z_half, nFilters_w_half);
		printf("Filters widths z=%d; w=%d\n", nFilters_z_half, nFilters_w_half);
		printf("Number of filters: %d\n", nFilters_total);
		printf("Halfwidth of the widest filter: %d\n", filter_halfwidth);
		printf("Useful part of the segment: %d\n", useful_part_size);
		printf("nSegments: %d\n", nSegments);
		printf("Number of z-planes: %d; Number of z-planes per chunk: %d;\n", nZPlanes, nZPlanes_per_chunk);
		printf("ZW chunks:\n");
		for(int f=0; f<(int) ZW_chunks.size(); f++){
			printf("    %d\n", ZW_chunks[f]);
		}
		printf("Number of chunks: %d; Remainder: %d;\n", nZPlanes_chunks, nZPlanes_remainder);
		printf("-------------------------------------------\n");
	}
};



#endif