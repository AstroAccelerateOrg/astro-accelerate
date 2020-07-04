#ifndef ASTRO_ACCELERATE_AA_JERK_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_JERK_STRATEGY_HPP

#include "aa_strategy.hpp"
#include "aa_jerk_plan.h"
#include "presto_funcs.hpp"
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
	size_t c_nSamples_time_dom;
	size_t c_nSamples_freq_dom;
	size_t c_nSamples_output_plane;
	size_t c_nDMs;
	
	// filter parameters
	float c_z_max_search_limit;
	float c_z_search_step;
	float c_w_max_search_limit;
	float c_w_search_step;
	int   c_interbinned_samples;
	int   c_high_precision;
	
	// MSD
	bool  c_MSD_outlier_rejection;
	float c_OR_sigma_cuttoff;
	
	// Candidate selection
	float c_CS_sigma_threshold;
	int   c_CS_algorithm;
	
	//-----------------> JERK strategy
	
	int c_conv_size; // size of the segment in overlap and save
	int c_nFilters_z_half;
	int c_nFilters_z;
	int c_nFilters_w_half;
	int c_nFilters_w;
	int c_nFilters_total;
	
	int c_filter_halfwidth;
	int c_useful_part_size;
	
	int    c_nSegments;
	size_t c_output_size_one_DM;
	size_t c_output_size_z_plane;
	size_t c_output_size_total;
	
	int c_nZPlanes;
	int c_nZPlanes_per_chunk;
	int c_nZPlanes_chunks;
	int c_nZPlanes_remainder;
	std::vector<int> c_ZW_chunks;
	size_t c_free_memory;
	size_t c_required_memory;
	size_t c_total_memory;
	size_t c_reserved_memory_for_candidate_selection;
	
	size_t c_filter_padded_size;
	size_t c_filter_padded_size_bytes;
	
	bool c_ready;
	
	bool calculate_memory_split(size_t available_free_memory){
		//---------> Memory testing and splitting
		c_nZPlanes = c_nFilters_w;
		size_t z_plane_size_bytes = c_output_size_z_plane*sizeof(float);
		available_free_memory = available_free_memory - c_reserved_memory_for_candidate_selection - filter_padded_size_bytes;
		c_nZPlanes_per_chunk = (int) (available_free_memory/z_plane_size_bytes);
		if(c_nZPlanes_per_chunk == 0) return(false);
		c_nZPlanes_chunks = (int) (nZPlanes/c_nZPlanes_per_chunk);
		c_nZPlanes_remainder = nZPlanes - c_nZPlanes_chunks*c_nZPlanes_per_chunk;
		
		for(int f=0; f<c_nZPlanes_chunks; f++){
			c_ZW_chunks.push_back(c_nZPlanes_per_chunk);
		}
		if(c_nZPlanes_remainder>0) {
			c_ZW_chunks.push_back(c_nZPlanes_remainder);
		}
		
		c_required_memory = (c_nZPlanes_chunks + c_nZPlanes_remainder)*c_z_plane_size_bytes + c_reserved_memory_for_candidate_selection;
		c_free_memory = available_free_memory - c_required_memory;
		return(true);
	}
	
public:
	aa_jerk_strategy(aa_jerk_plan plan, size_t available_free_memory){
		c_total_memory = available_free_memory;
		c_nSamples_time_dom   = plan.nTimesamples();
		c_nSamples_freq_dom   = (plan.nTimesamples()>>1) + 1; //because R2C FFT
		c_nDMs                = plan.nDMs();
		
		// Calculate number of filters
		// number of filters must also account for negative accelerations and w=z=0;
		c_nFilters_z_half   = plan.z_max_search_limit()/plan.z_search_step();
		c_nFilters_z        = c_nFilters_z_half + c_nFilters_z_half + 1; 
		c_nFilters_w_half   = plan.w_max_search_limit()/plan.w_search_step();
		c_nFilters_w        = c_nFilters_w_half + c_nFilters_w_half + 1;
		c_nFilters_total    = c_nFilters_z*c_nFilters_w;
		
		// recompute maximum z and w values based on step
		c_z_max_search_limit  = c_nFilters_z_half*plan.z_search_step();
		c_z_search_step       = plan.z_search_step();
		c_w_max_search_limit  = c_nFilters_w_half*plan.w_search_step();
		c_w_search_step       = plan.w_search_step();
		c_interbinned_samples = plan.number_of_interbinned_samples();
		c_high_precision      = plan.precision();
		
		// MSD
		c_MSD_outlier_rejection = plan.MSD_outlier_rejection();
		c_OR_sigma_cuttoff = plan.OR_sigma_cuttoff();
	
		// Candidate selection
		c_CS_sigma_threshold = plan.CS_sigma_threshold();
		c_CS_algorithm = plan.CS_algorithm();
		
		// Strategy
		c_conv_size = 2048;
		
		int presto_halfwidth = presto_w_resp_halfwidth(c_z_max_search_limit, c_w_max_search_limit, c_high_precision);
		c_filter_halfwidth = presto_halfwidth*c_interbinned_samples;
		c_useful_part_size = c_conv_size - 2*c_filter_halfwidth + 1;
		
		c_nSegments           = (c_nSamples_freq_dom + c_useful_part_size - 1)/c_useful_part_size;
		c_output_size_one_DM  = c_nSegments*c_useful_part_size;
		c_output_size_z_plane = c_nFilters_z*c_output_size_one_DM;
		c_output_size_total   = c_nFilters_total*c_output_size_one_DM;
		c_reserved_memory_for_candidate_selection = 2*c_output_size_z_plane*sizeof(float);
		
		c_filter_padded_size       = c_nFilters_total*c_conv_size;
		c_filter_padded_size_bytes = c_nFilters_total*c_conv_size*sizeof(float2); //*8 for complex float
		
		c_ready = calculate_memory_split(available_free_memory);
	}
	
	~aa_jerk_strategy(){
		c_ZW_chunks.clear();
	}
	
    bool ready() const {
      return c_ready;
    }
	
    bool setup() {
      return ready();
    }
	
	size_t nSamples_time_dom() return(c_nSamples_time_dom);
	size_t nSamples_freq_dom() return(c_nSamples_freq_dom);
	size_t nSamples_output_plane() return(c_nSamples_output_plane);
	size_t nDMs() return(c_nDMs);
	
	float z_max_search_limit() return(c_z_max_search_limit);
	float z_search_step() return(c_z_search_step);
	float w_max_search_limit() return(c_w_max_search_limit);
	float w_search_step() return(c_w_search_step);
	int   interbinned_samples() return(c_interbinned_samples);
	int   high_precision() return(c_high_precision);
	
	bool  MSD_outlier_rejection() return(c_MSD_outlier_rejection);
	float OR_sigma_cuttoff() return(c_OR_sigma_cuttoff);
	
	float CS_sigma_threshold() return(c_CS_sigma_threshold);
	int   CS_algorithm() return(c_CS_algorithm);
	
	int conv_size() return(c_conv_size);
	int nFilters_z_half() return(c_nFilters_z_half);
	int nFilters_z() return(c_nFilters_z);
	int nFilters_w_half() return(c_nFilters_w_half);
	int nFilters_w() return(c_nFilters_w);
	int nFilters_total() return(c_nFilters_total);
	
	int filter_halfwidth() return(c_filter_halfwidth);
	int useful_part_size() return(c_useful_part_size);
	
	int    nSegments() return(c_nSegments);
	size_t output_size_one_DM() return(c_output_size_one_DM);
	size_t output_size_z_plane() return(c_output_size_z_plane);
	size_t output_size_total() return(c_output_size_total);
	
	int nZPlanes() return(c_nZPlanes);
	int nZPlanes_per_chunk() return(c_nZPlanes_per_chunk);
	int nZPlanes_chunks() return(c_nZPlanes_chunks);
	int nZPlanes_remainder() return(c_nZPlanes_remainder);
	std::vector<int> ZW_chunks() return(c_ZW_chunks);
	size_t free_memory() return(c_free_memory);
	size_t required_memory() return(c_required_memory);
	size_t total_memory() return(c_total_memory);
	size_t reserved_memory_for_candidate_selection() return(c_reserved_memory_for_candidate_selection);
	
	size_t c_filter_padded_size() return(c_filter_padded_size);
	size_t c_filter_padded_size_bytes() return(c_filter_padded_size_bytes);
	
	void PrintStrategy(){
		printf("-------------------------------------------\n");
		printf("Input parameters:\n");
		printf("    Number of time samples before FFT: %zu\n", c_nSamples_time_dom);
		printf("    Number of time samples after FFT:  %zu\n", c_nSamples_freq_dom);
		printf("    Number of time samples in output:  %zu\n", c_nSamples_freq_dom);
		printf("    Number of DM trials:               %zu\n", c_nDMs);
		printf("Filter parameters:\n");
		printf("    Filter's halfwidth %d;\n", c_filter_halfwidth);
		printf("    z max:             %f;\n", c_z_max_search_limit);
		printf("    z step size:       %f;\n", c_z_search_step);
		printf("    w max:             %f;\n", c_w_max_search_limit);
		printf("    w step size:       %f;\n", c_w_search_step);
		printf("\n");
		printf("Interbinning: ");
		if(c_interbinned_samples==2) printf("Yes.\n"); else printf("No.\n");
		printf("High precision filters: ");
		if(c_high_precision==1) printf("Yes.\n"); else printf("No.\n");
		printf("-------------------------------------------\n");
		printf("\n");
		printf("-------------------------------------------\n");
		printf("Convolution size: %d\n", c_conv_size);
		printf("Half filters widths z=%d; w=%d\n", c_nFilters_z_half, c_nFilters_w_half);
		printf("Filters widths z=%d; w=%d\n", c_nFilters_z_half, c_nFilters_w_half);
		printf("Number of filters: %d\n", c_nFilters_total);
		printf("Halfwidth of the widest filter: %d\n", c_filter_halfwidth);
		printf("Useful part of the segment: %d\n", c_useful_part_size);
		printf("nSegments: %d\n", c_nSegments);
		printf("Number of z-planes: %d; Number of z-planes per chunk: %d;\n", c_nZPlanes, c_nZPlanes_per_chunk);
		printf("ZW chunks:\n");
		for(int f=0; f<(int) c_ZW_chunks.size(); f++){
			printf("    %d\n", c_ZW_chunks[f]);
		}
		printf("Number of chunks: %d; Remainder: %d;\n", c_nZPlanes_chunks, c_nZPlanes_remainder);
		printf("-------------------------------------------\n");
	}
};



#endif