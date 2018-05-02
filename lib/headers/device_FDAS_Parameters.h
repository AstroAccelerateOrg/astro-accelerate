#ifndef __ASTROACCELERATE_FDAS_PARAMETERS__
#define __ASTROACCELERATE_FDAS_PARAMETERS__

class FDAS_Parameters {
public:
	// Constants
	float sigma_cutoff;
	// Switches
	int candidate_algorithm;
	int enable_output_ffdot_plan;
	int enable_output_fdas_list;
	int enable_fdas_custom_fft;
	int enable_fdas_inbin;
	int enable_fdas_norm;
	
	FDAS_Parameters(){
		sigma_cutoff = 6.0f;
		candidate_algorithm = 0;
		enable_output_ffdot_plan = 0;
		enable_output_fdas_list = 0;
		enable_fdas_custom_fft = 0;
		enable_fdas_inbin = 0;
		enable_fdas_norm = 0;
	}
	
	void Set(float t_sigma_cutoff, int t_candidate_algorithm, int t_enable_output_ffdot_plan, int t_enable_output_fdas_list, int t_enable_fdas_custom_fft, int t_enable_fdas_inbin, int t_enable_fdas_norm){
		sigma_cutoff = t_sigma_cutoff;
		
		candidate_algorithm      = t_candidate_algorithm;
		enable_output_ffdot_plan = t_enable_output_ffdot_plan;
		enable_output_fdas_list  = t_enable_output_fdas_list;
		enable_fdas_custom_fft   = t_enable_fdas_custom_fft;
		enable_fdas_inbin        = t_enable_fdas_inbin;
		enable_fdas_norm         = t_enable_fdas_norm;
	}
};

#endif
