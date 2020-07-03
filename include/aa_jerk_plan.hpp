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
class aa_jerk_plan {
private:
	// input parameters
	size_t c_nTimesamples;
	size_t c_nDMs;
	
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
	
	// MSD
	bool  c_MSD_outlier_rejection;
	float c_OR_sigma_cuttoff;

public:
	
	aa_jerk_plan(){
		// input parameters
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
		c_OR_sigma_cuttoff = 0;
		
		// Candidate selection
		c_CS_sigma_threshold = 0;
		c_CS_algorithm = 0;
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
		c_nTimesamples = t_nTimesamples;
		c_nDMs         = t_nDMs;
		
		c_z_max_search_limit = t_z_max_search_limit;
		c_w_max_search_limit = t_w_max_search_limit;
		c_z_search_step      = t_z_search_step;
		c_w_search_step      = t_w_search_step;
		
		if(t_do_interbinning) c_interbinned_samples = 2; 
		else c_interbinned_samples = 1;
		
		if(t_do_high_precision) c_high_precision = 1;
		else c_high_precision = 0;
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
	
	size_t nTimesamples() return(c_nTimesamples);
	size_t nDMs() return(c_nDMs);
	float z_max_search_limit() return(c_z_max_search_limit);
	float z_search_step() return(c_z_search_step);
	float w_max_search_limit() return(c_w_max_search_limit);
	float w_search_step() return(c_w_search_step);
	int interbinned_samples() return(c_interbinned_samples);
	int high_precision() return(c_high_precision);
	
	float CS_sigma_threshold() return(c_CS_sigma_threshold);
	int CS_algorithm() return(c_CS_algorithm);
	
	bool MSD_outlier_rejection() return(c_MSD_outlier_rejection);
	float OR_sigma_cutoff() return(c_OR_sigma_cuttoff);
	
	void PrintPlan(){
		printf("Input parameters for FDAS/JERK search:\n");
		printf("    Number of time samples: %zu\n", c_nTimesamples);
		printf("    Number of DM trials:    %zu\n", c_nDMs);
		
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
	}
};

} // namespace astroaccelerate

#endif