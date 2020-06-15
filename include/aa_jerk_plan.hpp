#ifndef ASTRO_ACCELERATE_AA_JERK_PLAN_HPP
#define ASTRO_ACCELERATE_AA_JERK_PLAN_HPP


// NOTES:
// conv_size should really be determined by the class based on size of the filter and performance of the GPU
namespace astroaccelerate {

/** \class aa_jerk_plan aa_jerk_plan.hpp "include/aa_jerk_plan.hpp"
 * \brief Class which contains parameters of the jerk search.
 * \details An jerk search plan is where user define search parameters. It is required for creating jerk strategy.   
 * \author AstroAccelerate.
 * \date 2020-06-15
 */
class jerk_plan {
private:

public:
	// input parameters
	size_t nTimesamples;
	size_t nDMs;
	
	// filter parameters
	int filter_length;
	int filter_halfwidth;
	
	// filter parameters
	float z_max_search_limit;
	float z_search_step;
	float w_max_search_limit;
	float w_search_step;
	int interbinned_samples;
	int high_precision;
	
	// Candidate selection
	bool MSD_outlier_rejection;
	float OR_sigma_cuttoff;
	float CS_sigma_threshold;
	int CS_algorithm;
	
	aa_jerk_plan(){
		// input parameters
		nTimesamples = 0;
		nDMs = 0;
	
		// filter parameters
		filter_length = 0;
		filter_halfwidth = 0;
	
		// filter parameters
		z_max_search_limit = 0;
		z_search_step = 0;
		w_max_search_limit = 0;
		w_search_step = 0;
		interbinned_samples = 0;
		high_precision = 0;
	
		// Candidate selection
		MSD_outlier_rejection = false;
		OR_sigma_cuttoff = 0;
		CS_sigma_threshold = 0;
		CS_algorithm = 0;
	}
	
	aa_jerk_plan(
			size_t t_nTimesamples, 
			size_t t_nDMs,
			float t_z_max_search_limit, 
			float t_z_search_step, 
			float t_w_max_search_limit, 
			float t_w_search_step, 
			bool t_do_interbinning,
			bool t_do_high_precision){
		nTimesamples = t_nTimesamples;
		nDMs         = t_nDMs;
		
		z_max_search_limit = t_z_max_search_limit;
		w_max_search_limit = t_w_max_search_limit;
		z_search_step      = t_z_search_step;
		w_search_step      = t_w_search_step;
		
		if(t_do_interbinning) interbinned_samples = 2; 
		else interbinned_samples = 1;
		
		if(t_do_high_precision) high_precision = 1;
		else high_precision = 0;
	}
	
	void enable_MSD_outlier_rejection(){
		MSD_outlier_rejection = true;
	}
	
	void disable_MSD_outlier_rejection(){
		MSD_outlier_rejection = false;
	}
	
	void set_outlier_rejection_sigma_cutoff(float sigma) {
		OR_sigma_cuttoff = sigma;
	}
	
	void set_candidate_selection_sigma_threshold(float sigma) {
		CS_sigma_threshold = sigma;
	}
	
	void PrintPlan(){
		printf("Input parameters for FDAS/JERK search:\n");
		printf("    Number of time samples: %zu\n", nTimesamples);
		printf("    Number of DM trials:    %zu\n", nDMs);
		
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
	}
};

} // namespace astroaccelerate

#endif