#ifndef ASTRO_ACCELERATE_AA_JERK_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_JERK_STRATEGY_HPP

#include <stdio.h>

#include "aa_params.hpp"
#include "aa_strategy.hpp"
#include "aa_jerk_plan.hpp"
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
	int64_t c_nSamples_time_dom; // length of the time-domain input signal to FFT
	int64_t c_nSamples_freq_dom; // length of the freq-domain output signal from FFT
	int64_t c_nSamples_output_plane; // Does nothing ...
	int64_t c_nTimesamples; // number of timesamples of DM trial
	int64_t c_nDMs; // Number of DM trials
	
	// filter parameters
	float   c_z_max_search_limit;
	float   c_z_search_step;
	float   c_w_max_search_limit;
	float   c_w_search_step;
	int64_t c_interbinned_samples;
	int64_t c_high_precision;
	
	// MSD
	bool    c_MSD_outlier_rejection;
	float   c_OR_sigma_cutoff;
	
	// Candidate selection
	float   c_CS_sigma_threshold;
	int64_t c_CS_algorithm;
	int64_t c_nHarmonics;
	
	// other flags
	bool c_always_choose_next_power_of_2;
	bool c_spectrum_whitening;
	
	//-----------------> JERK strategy
	
	int64_t c_conv_size; // size of the segment in overlap and save
	int64_t c_nFilters_z_half;
	int64_t c_nFilters_z;
	int64_t c_nFilters_w_half;
	int64_t c_nFilters_w;
	int64_t c_nFilters_total;
	
	int64_t c_filter_halfwidth;
	int64_t c_useful_part_size;
	
	int c_nSegments;
	int64_t c_output_size_one_DM;
	int64_t c_output_size_z_plane;
	int64_t c_output_size_total;
	
	int64_t c_nZPlanes;
	int64_t c_nZPlanes_per_chunk;
	int64_t c_nZPlanes_chunks;
	int64_t c_nZPlanes_remainder;
	std::vector<int64_t> c_ZW_chunks;
	int64_t c_free_memory;
	int64_t c_required_memory;
	int64_t c_available_memory;
	int64_t c_reserved_memory_for_candidate_selection;
	
	int64_t c_filter_padded_size;
	int64_t c_filter_padded_size_bytes;
	
	bool c_ready;
	
	int64_t next_power_2(int64_t x) {
		x--;
		x |= x>>1;
		x |= x>>2;
		x |= x>>4;
		x |= x>>8;
		x |= x>>16;
		x |= x>>32;
		x++;
		return x;
	}
	
	bool calculate_memory_split(int64_t available_free_memory){
		//---------> Memory testing and splitting
		c_nZPlanes = c_nFilters_w;
		int64_t z_plane_size_bytes = c_output_size_z_plane*sizeof(float);
		int64_t MSD_workarea_size = ((int) ((c_nFilters_z*c_nTimesamples)/PD_NTHREADS))*MSD_PARTIAL_SIZE*sizeof(float);
		available_free_memory = available_free_memory - 2*c_reserved_memory_for_candidate_selection - c_filter_padded_size_bytes - c_nSamples_time_dom*5*sizeof(float) - MSD_workarea_size;
		c_nZPlanes_per_chunk = (int) (available_free_memory/z_plane_size_bytes);
		if(c_nZPlanes_per_chunk == 0) return(false);
		c_nZPlanes_chunks = (int) (c_nZPlanes/c_nZPlanes_per_chunk);
		c_nZPlanes_remainder = c_nZPlanes - c_nZPlanes_chunks*c_nZPlanes_per_chunk;
		
		for(int f=0; f<c_nZPlanes_chunks; f++){
			c_ZW_chunks.push_back(c_nZPlanes_per_chunk);
		}
		if(c_nZPlanes_remainder>0) {
			c_ZW_chunks.push_back(c_nZPlanes_remainder);
		}
		
		if(c_ZW_chunks.size()>0) {
			c_required_memory = c_ZW_chunks[0]*z_plane_size_bytes;
		}
		else return (false);
		
		c_free_memory = available_free_memory - c_required_memory;
		return(true);
	}
	
public:
	aa_jerk_strategy() {
		c_nSamples_time_dom = 0;
		c_nSamples_freq_dom = 0;
		c_nSamples_output_plane = 0;
		c_nTimesamples = 0;
		c_nDMs = 0;
		
		// filter parameters
		c_z_max_search_limit = 0;
		c_z_search_step = 0;
		c_w_max_search_limit = 0;
		c_w_search_step = 0;
		c_interbinned_samples = 0;
		c_high_precision = 0;
		
		// MSD
		c_MSD_outlier_rejection = false;
		c_OR_sigma_cutoff = 0;
		
		// Candidate selection
		c_CS_sigma_threshold = 0;
		c_CS_algorithm = 1;
		c_nHarmonics = 0;
		
		// other flags
		c_always_choose_next_power_of_2 = false;
		c_spectrum_whitening = false;
		
		c_conv_size = 2048;
		c_nFilters_z_half = 0;
		c_nFilters_z = 0;
		c_nFilters_w_half = 0;
		c_nFilters_w = 0;
		c_nFilters_total = 0;
		
		c_filter_halfwidth = 0;
		c_useful_part_size = 0;
		
		c_nSegments = 0;
		c_output_size_one_DM = 0;
		c_output_size_z_plane = 0;
		c_output_size_total = 0;
		
		c_nZPlanes = 0;
		c_nZPlanes_per_chunk = 0;
		c_nZPlanes_chunks = 0;
		c_nZPlanes_remainder = 0;
		c_free_memory = 0;
		c_required_memory = 0;
		c_available_memory = 0;
		c_reserved_memory_for_candidate_selection = 0;
		
		c_filter_padded_size = 0;
		c_filter_padded_size_bytes = 0;
		
		c_ready = false;
	}
	
	aa_jerk_strategy(aa_jerk_plan &plan){
		c_available_memory = plan.available_memory();
		c_always_choose_next_power_of_2 = plan.always_choose_next_power_of_2();
		c_spectrum_whitening = plan.spectrum_whitening();
		
		c_nTimesamples      = plan.nTimesamples();
		c_nSamples_time_dom = next_power_2(c_nTimesamples);
		if(!c_always_choose_next_power_of_2) 
			c_nSamples_time_dom = (c_nSamples_time_dom>>1);
		c_nSamples_freq_dom = (c_nSamples_time_dom>>1) + 1; //because R2C FFT
		c_nDMs              = plan.nDMs();
		
		// Calculate number of filters
		// number of filters must also account for negative accelerations and w=z=0;
		if(plan.z_search_step()>0) 
			c_nFilters_z_half = plan.z_max_search_limit()/plan.z_search_step();
		else c_nFilters_z_half = 0;
		c_nFilters_z        = c_nFilters_z_half + c_nFilters_z_half + 1; 
		
		if(plan.w_search_step()>0) 
			c_nFilters_w_half = plan.w_max_search_limit()/plan.w_search_step();
		else c_nFilters_w_half = 0;
		c_nFilters_w      = c_nFilters_w_half + c_nFilters_w_half + 1;
		
		c_nFilters_total  = c_nFilters_z*c_nFilters_w;
		
		// recompute maximum z and w values based on step
		c_z_max_search_limit  = c_nFilters_z_half*plan.z_search_step();
		c_z_search_step       = plan.z_search_step();
		c_w_max_search_limit  = c_nFilters_w_half*plan.w_search_step();
		c_w_search_step       = plan.w_search_step();
		c_interbinned_samples = plan.interbinned_samples();
		c_high_precision      = plan.high_precision();
		
		// MSD
		c_MSD_outlier_rejection = plan.MSD_outlier_rejection();
		c_OR_sigma_cutoff = plan.OR_sigma_cutoff();
		c_nHarmonics = plan.nHarmonics();
	
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
		
		c_ready = calculate_memory_split(c_available_memory);
	}
	
	bool recalculate(int64_t nTimesamples, int64_t nDMs){
		c_ready = false;
		c_ZW_chunks.clear();
		c_nTimesamples      = nTimesamples;
		c_nSamples_time_dom = next_power_2(c_nTimesamples);
		if(!c_always_choose_next_power_of_2) 
			c_nSamples_time_dom = (c_nSamples_time_dom>>1);
		c_nSamples_freq_dom = (c_nSamples_time_dom>>1) + 1; //because R2C FFT
		c_nDMs              = nDMs;
		
		c_nSegments           = (c_nSamples_freq_dom + c_useful_part_size - 1)/c_useful_part_size;c_output_size_one_DM  = c_nSegments*c_useful_part_size;
		c_output_size_z_plane = c_nFilters_z*c_output_size_one_DM;
		c_output_size_total   = c_nFilters_total*c_output_size_one_DM;
		c_reserved_memory_for_candidate_selection = 2*c_output_size_z_plane*sizeof(float);
		
		c_ready = calculate_memory_split(c_available_memory);
		return(c_ready);
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
	
	/** \returns The name of this mdoule. */
    std::string name() const {
      return "jerk_strategy";
    }
	
	int64_t nSamples_time_dom() {return(c_nSamples_time_dom);}
	int64_t nSamples_freq_dom() {return(c_nSamples_freq_dom);}
	int64_t nSamples_output_plane() {return(c_nSamples_output_plane);}
	int64_t nDMs() {return(c_nDMs);}
	int64_t nTimesamples() {return(c_nTimesamples);}
	
	float z_max_search_limit() {return(c_z_max_search_limit);}
	float z_search_step() {return(c_z_search_step);}
	float w_max_search_limit() {return(c_w_max_search_limit);}
	float w_search_step() {return(c_w_search_step);}
	int64_t interbinned_samples() {return(c_interbinned_samples);}
	int64_t high_precision() {return(c_high_precision);}
	
	bool  MSD_outlier_rejection() {return(c_MSD_outlier_rejection);}
	float OR_sigma_cutoff() {return(c_OR_sigma_cutoff);}
	
	float CS_sigma_threshold() {return(c_CS_sigma_threshold);}
	int64_t CS_algorithm() {return(c_CS_algorithm);}
	int64_t nHarmonics() {return(c_nHarmonics);}
	
	bool always_choose_next_power_of_2() {return(c_always_choose_next_power_of_2);}
	bool spectrum_whitening() {return(c_spectrum_whitening);}
	
	int64_t conv_size() {return(c_conv_size);}
	int64_t nFilters_z_half() {return(c_nFilters_z_half);}
	int64_t nFilters_z() {return(c_nFilters_z);}
	int64_t nFilters_w_half() {return(c_nFilters_w_half);}
	int64_t nFilters_w() {return(c_nFilters_w);}
	int64_t nFilters_total() {return(c_nFilters_total);}
	
	int64_t filter_halfwidth() {return(c_filter_halfwidth);}
	int64_t useful_part_size() {return(c_useful_part_size);}
	
	int64_t nSegments() {return(c_nSegments);}
	int64_t output_size_one_DM() {return(c_output_size_one_DM);}
	int64_t output_size_z_plane() {return(c_output_size_z_plane);}
	int64_t output_size_total() {return(c_output_size_total);}
	
	int64_t nZPlanes() {return(c_nZPlanes);}
	int64_t nZPlanes_per_chunk() {return(c_nZPlanes_per_chunk);}
	int64_t nZPlanes_chunks() {return(c_nZPlanes_chunks);}
	int64_t nZPlanes_remainder() {return(c_nZPlanes_remainder);}
	std::vector<int64_t> ZW_chunks() {return(c_ZW_chunks);}
	int64_t free_memory() {return(c_free_memory);}
	int64_t required_memory() {return(c_required_memory);}
	int64_t available_memory() {return(c_available_memory);}
	int64_t reserved_memory_for_candidate_selection() {return(c_reserved_memory_for_candidate_selection);}
	
	int64_t filter_padded_size() {return(c_filter_padded_size);}
	int64_t filter_padded_size_bytes() {return(c_filter_padded_size_bytes);}
	
	static bool print_info(aa_jerk_strategy &strategy){
		printf("-------------------------------------------\n");
		printf("Input parameters:\n");
		printf("    Number of time samples:            %ld\n", strategy.nTimesamples());
		printf("    Number of time samples (padded):   %ld\n", strategy.nSamples_time_dom());
		printf("    Number of time samples after FFT:  %ld\n", strategy.nSamples_freq_dom());
		printf("    Number of time samples in output:  %ld\n", strategy.nSamples_freq_dom());
		printf("    Number of DM trials:               %ld\n", strategy.nDMs());
		printf("Filter parameters:\n");
		printf("    Filter's halfwidth %ld;\n", strategy.filter_halfwidth());
		printf("    z max:             %f;\n", strategy.z_max_search_limit());
		printf("    z step size:       %f;\n", strategy.z_search_step());
		printf("    w max:             %f;\n", strategy.w_max_search_limit());
		printf("    w step size:       %f;\n", strategy.w_search_step());
		printf("\n");
		printf("Interbinning: ");
		if(strategy.interbinned_samples()==2) printf("Yes.\n"); else printf("No.\n");
		printf("High precision filters: ");
		if(strategy.high_precision()==1) printf("Yes.\n"); else printf("No.\n");
		printf("\n");
		printf("Candidate selection:\n");
		printf("    CS sigma thresh:   %f;\n", strategy.CS_sigma_threshold());
		printf("    Number harmonics:  %ld;\n", strategy.nHarmonics());
		printf("    Algoriths:         ");
		if(strategy.CS_algorithm()==0) printf("threshold\n");
		else if(strategy.CS_algorithm()==1) printf("peak finding\n");
		printf("Mean and standard deviation:\n");
		printf("    Outlier rejection: ");
		if(strategy.MSD_outlier_rejection()==1) printf("Yes.\n"); else printf("No.\n");
		printf("    Sigma cutoff:      %f;\n", strategy.OR_sigma_cutoff());
		printf("Other flags:\n");
		printf("    Next power of two:      ");
		if(strategy.always_choose_next_power_of_2()) printf("Yes.\n"); else printf("No.\n");
		printf("    Spectrum whitening:     ");
		if(strategy.spectrum_whitening()) printf("Yes.\n"); else printf("No.\n");
		printf("-------------------------------------------\n");
		printf("\n");
		printf("-------------------------------------------\n");
		printf("Convolution size: %ld\n", strategy.conv_size());
		printf("Number of filters in positive half z=%ld; w=%ld\n", strategy.nFilters_z_half(), strategy.nFilters_w_half());
		printf("Number of filters in z=%ld; w=%ld\n", strategy.nFilters_z(), strategy.nFilters_w());
		printf("Total number of filters: %ld\n", strategy.nFilters_total());
		printf("Halfwidth of the widest filter: %ld\n", strategy.filter_halfwidth());
		printf("Useful part of the segment: %ld\n", strategy.useful_part_size());
		printf("nSegments: %ld\n", strategy.nSegments());
		printf("Number of z-planes: %ld; Number of z-planes per chunk: %ld;\n", strategy.nZPlanes(), strategy.nZPlanes_per_chunk());
		printf("ZW chunks:\n");
		std::vector<int64_t> ZW_chunks = strategy.ZW_chunks();
		for(size_t f=0; f<ZW_chunks.size(); f++){
			printf("    %ld\n", ZW_chunks[f]);
		}
		printf("Number of chunks: %ld; Remainder: %ld;\n", strategy.nZPlanes_chunks(), strategy.nZPlanes_remainder());
		printf("-------------------------------------------\n");
		return true;
	}
};

} // name space

#endif