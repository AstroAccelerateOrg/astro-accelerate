#ifndef ASTRO_ACCELERATE_AA_JERK_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_JERK_STRATEGY_HPP

#include <stdio.h>

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
	unsigned long int c_nSamples_time_dom;
	unsigned long int c_nSamples_freq_dom;
	unsigned long int c_nSamples_output_plane;
	unsigned long int c_nTimesamples;
	unsigned long int c_nDMs;
	
	// filter parameters
	float c_z_max_search_limit;
	float c_z_search_step;
	float c_w_max_search_limit;
	float c_w_search_step;
	int   c_interbinned_samples;
	int   c_high_precision;
	
	// MSD
	bool  c_MSD_outlier_rejection;
	float c_OR_sigma_cutoff;
	
	// Candidate selection
	float c_CS_sigma_threshold;
	int   c_CS_algorithm;
	int   c_nHarmonics;
	
	// other flags
	bool c_always_choose_next_power_of_2;
	bool c_spectrum_whitening;
	
	//-----------------> JERK strategy
	
	int c_conv_size; // size of the segment in overlap and save
	int c_nFilters_z_half;
	int c_nFilters_z;
	int c_nFilters_w_half;
	int c_nFilters_w;
	int c_nFilters_total;
	
	int c_filter_halfwidth;
	int c_useful_part_size;
	
	int c_nSegments;
	unsigned long int c_output_size_one_DM;
	unsigned long int c_output_size_z_plane;
	unsigned long int c_output_size_total;
	
	int c_nZPlanes;
	int c_nZPlanes_per_chunk;
	int c_nZPlanes_chunks;
	int c_nZPlanes_remainder;
	std::vector<int> c_ZW_chunks;
	unsigned long int c_free_memory;
	unsigned long int c_required_memory;
	unsigned long int c_total_memory;
	unsigned long int c_reserved_memory_for_candidate_selection;
	
	unsigned long int c_filter_padded_size;
	unsigned long int c_filter_padded_size_bytes;
	
	bool c_ready;
	
	unsigned long int next_power_2(unsigned long int x) {
		/*
		int exp = 0;
		unsigned long int temp = x;
		while(temp!=0){
			temp = temp>>1;
			exp++;
		}
		printf("exp: %d; temp: %zu;\n", exp, temp);
		temp = 1<<exp;
		printf("temp: %zu\n", temp);
		return(temp);
		*/
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
	
	bool calculate_memory_split(unsigned long int available_free_memory){
		//---------> Memory testing and splitting
		c_nZPlanes = c_nFilters_w;
		unsigned long int z_plane_size_bytes = c_output_size_z_plane*sizeof(float);
		available_free_memory = available_free_memory - c_reserved_memory_for_candidate_selection - c_filter_padded_size_bytes;
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
			c_required_memory = c_ZW_chunks[0]*z_plane_size_bytes + c_reserved_memory_for_candidate_selection;
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
		c_total_memory = 0;
		c_reserved_memory_for_candidate_selection = 0;
		
		c_filter_padded_size = 0;
		c_filter_padded_size_bytes = 0;
		
		c_ready = false;
	}
	
	aa_jerk_strategy(aa_jerk_plan &plan){
		c_total_memory = plan.available_memory();
		c_always_choose_next_power_of_2 = plan.always_choose_next_power_of_2();
		c_spectrum_whitening = plan.spectrum_whitening();
		
		c_nTimesamples      = plan.nTimesamples();
		c_nSamples_time_dom = next_power_2(c_nTimesamples);
		if( (c_nTimesamples/c_nSamples_time_dom)<0.65 && !c_always_choose_next_power_of_2) c_nSamples_time_dom = (c_nSamples_time_dom>>1);
		c_nSamples_freq_dom = (c_nSamples_time_dom>>1) + 1; //because R2C FFT
		c_nDMs              = plan.nDMs();
		
		// Calculate number of filters
		// number of filters must also account for negative accelerations and w=z=0;
		if(plan.z_search_step()>0) c_nFilters_z_half = plan.z_max_search_limit()/plan.z_search_step();
		else c_nFilters_z_half = 0;
		c_nFilters_z        = c_nFilters_z_half + c_nFilters_z_half + 1; 
		
		if(plan.w_search_step()>0) c_nFilters_w_half = plan.w_max_search_limit()/plan.w_search_step();
		c_nFilters_w_half = 0;
		c_nFilters_w        = c_nFilters_w_half + c_nFilters_w_half + 1;
		
		c_nFilters_total    = c_nFilters_z*c_nFilters_w;
		
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
		
		c_ready = calculate_memory_split(c_total_memory);
	}
	
	bool recalculate(unsigned long int nTimesamples, unsigned long int nDMs){
		c_ready = false;
		c_ZW_chunks.clear();
		c_nTimesamples      = nTimesamples;
		c_nSamples_time_dom = next_power_2(c_nTimesamples);
		if( (c_nTimesamples/c_nSamples_time_dom)<0.65 && !c_always_choose_next_power_of_2) c_nSamples_time_dom = (c_nSamples_time_dom>>1);
		c_nSamples_freq_dom = (c_nSamples_time_dom>>1) + 1; //because R2C FFT
		c_nDMs              = nDMs;
		
		c_nSegments           = (c_nSamples_freq_dom + c_useful_part_size - 1)/c_useful_part_size;c_output_size_one_DM  = c_nSegments*c_useful_part_size;
		c_output_size_z_plane = c_nFilters_z*c_output_size_one_DM;
		c_output_size_total   = c_nFilters_total*c_output_size_one_DM;
		c_reserved_memory_for_candidate_selection = 2*c_output_size_z_plane*sizeof(float);
		
		c_ready = calculate_memory_split(c_total_memory);
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
	
	unsigned long int nSamples_time_dom() {return(c_nSamples_time_dom);}
	unsigned long int nSamples_freq_dom() {return(c_nSamples_freq_dom);}
	unsigned long int nSamples_output_plane() {return(c_nSamples_output_plane);}
	unsigned long int nDMs() {return(c_nDMs);}
	unsigned long int nTimesamples() {return(c_nTimesamples);}
	
	float z_max_search_limit() {return(c_z_max_search_limit);}
	float z_search_step() {return(c_z_search_step);}
	float w_max_search_limit() {return(c_w_max_search_limit);}
	float w_search_step() {return(c_w_search_step);}
	int   interbinned_samples() {return(c_interbinned_samples);}
	int   high_precision() {return(c_high_precision);}
	
	bool  MSD_outlier_rejection() {return(c_MSD_outlier_rejection);}
	float OR_sigma_cutoff() {return(c_OR_sigma_cutoff);}
	
	float CS_sigma_threshold() {return(c_CS_sigma_threshold);}
	int   CS_algorithm() {return(c_CS_algorithm);}
	int   nHarmonics() {return(c_nHarmonics);}
	
	bool always_choose_next_power_of_2() {return(c_always_choose_next_power_of_2);}
	bool spectrum_whitening() {return(c_spectrum_whitening);}
	
	int conv_size() {return(c_conv_size);}
	int nFilters_z_half() {return(c_nFilters_z_half);}
	int nFilters_z() {return(c_nFilters_z);}
	int nFilters_w_half() {return(c_nFilters_w_half);}
	int nFilters_w() {return(c_nFilters_w);}
	int nFilters_total() {return(c_nFilters_total);}
	
	int filter_halfwidth() {return(c_filter_halfwidth);}
	int useful_part_size() {return(c_useful_part_size);}
	
	int nSegments() {return(c_nSegments);}
	unsigned long int output_size_one_DM() {return(c_output_size_one_DM);}
	unsigned long int output_size_z_plane() {return(c_output_size_z_plane);}
	unsigned long int output_size_total() {return(c_output_size_total);}
	
	int nZPlanes() {return(c_nZPlanes);}
	int nZPlanes_per_chunk() {return(c_nZPlanes_per_chunk);}
	int nZPlanes_chunks() {return(c_nZPlanes_chunks);}
	int nZPlanes_remainder() {return(c_nZPlanes_remainder);}
	std::vector<int> ZW_chunks() {return(c_ZW_chunks);}
	unsigned long int free_memory() {return(c_free_memory);}
	unsigned long int required_memory() {return(c_required_memory);}
	unsigned long int total_memory() {return(c_total_memory);}
	unsigned long int reserved_memory_for_candidate_selection() {return(c_reserved_memory_for_candidate_selection);}
	
	unsigned long int filter_padded_size() {return(c_filter_padded_size);}
	unsigned long int filter_padded_size_bytes() {return(c_filter_padded_size_bytes);}
	
	static bool print_info(aa_jerk_strategy &strategy){
		printf("-------------------------------------------\n");
		printf("Input parameters:\n");
		printf("    Number of time samples:            %zu\n", strategy.nTimesamples());
		printf("    Number of time samples (padded):   %zu\n", strategy.nSamples_time_dom());
		printf("    Number of time samples after FFT:  %zu\n", strategy.nSamples_freq_dom());
		printf("    Number of time samples in output:  %zu\n", strategy.nSamples_freq_dom());
		printf("    Number of DM trials:               %zu\n", strategy.nDMs());
		printf("Filter parameters:\n");
		printf("    Filter's halfwidth %d;\n", strategy.filter_halfwidth());
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
		printf("    Number harmonics:  %d;\n", strategy.nHarmonics());
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
		printf("Convolution size: %d\n", strategy.conv_size());
		printf("Number of filters in positive half z=%d; w=%d\n", strategy.nFilters_z_half(), strategy.nFilters_w_half());
		printf("Number of filters in z=%d; w=%d\n", strategy.nFilters_z(), strategy.nFilters_w());
		printf("Number of filters: %d\n", strategy.nFilters_total());
		printf("Halfwidth of the widest filter: %d\n", strategy.filter_halfwidth());
		printf("Useful part of the segment: %d\n", strategy.useful_part_size());
		printf("nSegments: %d\n", strategy.nSegments());
		printf("Number of z-planes: %d; Number of z-planes per chunk: %d;\n", strategy.nZPlanes(), strategy.nZPlanes_per_chunk());
		printf("ZW chunks:\n");
		std::vector<int> ZW_chunks = strategy.ZW_chunks();
		for(int f=0; f<(int) ZW_chunks.size(); f++){
			printf("    %d\n", ZW_chunks[f]);
		}
		printf("Number of chunks: %d; Remainder: %d;\n", strategy.nZPlanes_chunks(), strategy.nZPlanes_remainder());
		printf("-------------------------------------------\n");
		return true;
	}
};

} // name space

#endif