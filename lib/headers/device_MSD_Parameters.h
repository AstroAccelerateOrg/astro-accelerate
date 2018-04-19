#ifndef __ASTROACCELERATE_MSD_PARAMETERS__
#define __ASTROACCELERATE_MSD_PARAMETERS__

class MSD_Parameters {
public:
	// Constants
	float OR_sigma_multiplier;
	// Switches
	int enable_outlier_rejection;
	
	MSD_Parameters(float t_OR_sigma_multiplier, int t_enable_outlier_rejection){
		OR_sigma_multiplier      = t_OR_sigma_multiplier;
		enable_outlier_rejection = t_enable_outlier_rejection;
	}
};

#endif
