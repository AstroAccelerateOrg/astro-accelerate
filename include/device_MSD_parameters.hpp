#ifndef ASTRO_ACCELERATE_MSD_PARAMETERS_HPP
#define ASTRO_ACCELERATE_MSD_PARAMETERS_HPP

class MSD_Parameters {
public:
	// Constants
	float OR_sigma_multiplier;
	// Switches
	int enable_outlier_rejection;
	
	MSD_Parameters(){
		enable_outlier_rejection = 0;
		OR_sigma_multiplier = 4.0f;
	}
	
	void set(float t_OR_sigma_multiplier, int t_enable_outlier_rejection){
		OR_sigma_multiplier      = t_OR_sigma_multiplier;
		enable_outlier_rejection = t_enable_outlier_rejection;
	}
	
	void debug(){
		printf("MSD parameters\n");
		printf("enable_outlier_rejection = %d;\n", enable_outlier_rejection);
		printf("OR_sigma_multiplier = %f;\n", OR_sigma_multiplier);
	}
};

#endif