#ifndef ASTRO_ACCELERATE_AA_JERK_PLAN_HPP
#define ASTRO_ACCELERATE_AA_JERK_PLAN_HPP

#include <stdio.h>

// NOTES:
// conv_size should really be determined by the class based on size of the filter and performance of the GPU
namespace astroaccelerate {

/** \class aa_jerk_plan aa_jerk_plan.hpp "include/aa_jerk_plan.hpp"
 * \brief Class which contains parameters of the jerk search.
 * \details An jerk search plan is where user define search parameters. It is required for creating jerk strategy.   
 * \author AstroAccelerate.
 * \date 2020-06-15
 */
class aa_jerk_plan {
private:
	// input parameters
	unsigned long int c_nTimesamples;
	unsigned long int c_nDMs;
	unsigned long int c_available_memory;
	
	// filter parameters
	float c_z_max_search_limit;
	float c_z_search_step;
	float c_w_max_search_limit;
	float c_w_search_step;
	int   c_interbinned_samples;
	int   c_high_precision;

	// Candidate selection
	float c_CS_sigma_threshold;
	int   c_CS_algorithm;
	int   c_nHarmonics;
	
	// MSD
	bool  c_MSD_outlier_rejection;
	float c_OR_sigma_cuttoff;
	
	// other flags
	bool c_always_choose_next_power_of_2;
	bool c_spectrum_whitening;
	

public:
	
	aa_jerk_plan(){
		// input parameters
		c_nTimesamples = 0;
		c_nDMs = 0;
		c_available_memory = 0;
	
		// filter parameters
		c_z_max_search_limit = 0;
		c_z_search_step = 0;
		c_w_max_search_limit = 0;
		c_w_search_step = 0;
		c_interbinned_samples = 0;
		c_high_precision = 0;
	
		// MSD
		c_MSD_outlier_rejection = false;
		c_OR_sigma_cuttoff = 0;
		
		// Candidate selection
		c_CS_sigma_threshold = 0;
		c_CS_algorithm = 0;
		c_nHarmonics = 8;
		
		// other flags
		c_always_choose_next_power_of_2 = false;
		c_spectrum_whitening = false;
	}
	
	aa_jerk_plan(
			unsigned long int nTimesamples, 
			unsigned long int nDMs,
			unsigned long int available_memory,
			float z_max_search_limit, 
			float z_search_step, 
			float w_max_search_limit, 
			float w_search_step, 
			bool do_interbinning,
			bool do_high_precision,
			float candidate_selection_sigma_threshold,
			int nHarmonics,
			bool do_outlier_rejection,
			float outlier_rejection_sigma_cutoff,
			bool always_choose_next_power_of_2,
			bool spectrum_whitening){
		c_nTimesamples = nTimesamples;
		c_nDMs         = nDMs;
		c_available_memory = available_memory;
		
		c_z_max_search_limit = z_max_search_limit;
		c_w_max_search_limit = w_max_search_limit;
		c_z_search_step      = z_search_step;
		c_w_search_step      = w_search_step;
		if(do_interbinning) c_interbinned_samples = 2; 
		else c_interbinned_samples = 1;
		if(do_high_precision) c_high_precision = 1;
		else c_high_precision = 0;
		
		c_CS_sigma_threshold = candidate_selection_sigma_threshold;
		c_CS_algorithm = 0;
		c_nHarmonics = nHarmonics;
		
		c_MSD_outlier_rejection = do_outlier_rejection;
		c_OR_sigma_cuttoff = outlier_rejection_sigma_cutoff;
		
		c_always_choose_next_power_of_2 = always_choose_next_power_of_2;
		c_spectrum_whitening = spectrum_whitening;
	}
	
	void change_nTimesamples(unsigned long int nTimesamples){
		c_nTimesamples = nTimesamples;
	}
	
	void change_nDMs(unsigned long int nDMs){
		c_nDMs = nDMs;
	}
	
	void enable_MSD_outlier_rejection(){
		c_MSD_outlier_rejection = true;
	}
	
	void disable_MSD_outlier_rejection(){
		c_MSD_outlier_rejection = false;
	}
	
	void enable_high_precision_filters(){
		c_high_precision = 1;
	}
	
	void disable_high_precision_filters(){
		c_high_precision = 0;
	}
	
	void enable_interbinning(){
		c_interbinned_samples = 2;
	}
	
	void disable_interbinning(){
		c_interbinned_samples = 1;
	}
	
	void set_outlier_rejection_sigma_cutoff(float sigma) {
		c_OR_sigma_cuttoff = sigma;
	}
	
	void set_candidate_selection_sigma_threshold(float sigma) {
		c_CS_sigma_threshold = sigma;
	}
	
	unsigned long int nTimesamples() {return(c_nTimesamples);}
	unsigned long int nDMs() {return(c_nDMs);}
	unsigned long int available_memory() {return(c_available_memory);}
	float z_max_search_limit(){return(c_z_max_search_limit);}
	float z_search_step() {return(c_z_search_step);}
	float w_max_search_limit() {return(c_w_max_search_limit);}
	float w_search_step() {return(c_w_search_step);}
	int interbinned_samples() {return(c_interbinned_samples);}
	int high_precision() {return(c_high_precision);}
	
	float CS_sigma_threshold() {return(c_CS_sigma_threshold);}
	int CS_algorithm() {return(c_CS_algorithm);}
	int nHarmonics() {return(c_nHarmonics);}
	
	bool MSD_outlier_rejection() {return(c_MSD_outlier_rejection);}
	float OR_sigma_cutoff() {return(c_OR_sigma_cuttoff);}
	
	bool always_choose_next_power_of_2() {return(c_always_choose_next_power_of_2);}
	bool spectrum_whitening() {return(c_spectrum_whitening);}
	
	void PrintPlan(){
		printf("Input parameters for FDAS/JERK search:\n");
		printf("    Number of time samples: %zu\n", c_nTimesamples);
		printf("    Number of DM trials:    %zu\n", c_nDMs);
		printf("    Available memory:       %zu\n", c_available_memory);
		
		printf("Filter parameters:\n");
		printf("    z max:             %f;\n", c_z_max_search_limit);
		printf("    z step size:       %f;\n", c_z_search_step);
		printf("    w max:             %f;\n", c_w_max_search_limit);
		printf("    w step size:       %f;\n", c_w_search_step);
		printf("\n");
		printf("Interbinning: ");
		if(c_interbinned_samples==2) printf("Yes.\n"); else printf("No.\n");
		printf("High precision filters: ");
		if(c_high_precision==1) printf("Yes.\n"); else printf("No.\n");
		
		printf("    CS sigma thresh:   %f;\n", c_CS_sigma_threshold);
		printf("    Number harmonics:  %d;\n", c_nHarmonics);
		printf("    Algoriths:         ");
		if(c_CS_algorithm==0) printf("threshold\n");
		else if(c_CS_algorithm==1) printf("peak finding\n");
		
		printf("Outlier rejection     : ");
		if(c_MSD_outlier_rejection==1) printf("Yes.\n"); else printf("No.\n");
		printf("    Sigma cutoff:      %f;\n", c_OR_sigma_cuttoff);
		
		printf("Next power of two:      ");
		if(c_always_choose_next_power_of_2) printf("Yes.\n"); else printf("No.\n");
		printf("Spectrum whitening:     ");
		if(c_spectrum_whitening) printf("Yes.\n"); else printf("No.\n");
	}
};

} // namespace astroaccelerate

#endif