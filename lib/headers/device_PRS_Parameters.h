#ifndef __ASTROACCELERATE_PRS_PARAMETERS__
#define __ASTROACCELERATE_PRS_PARAMETERS__

class PRS_Parameters {
public:
	// Constants
	float sigma_cutoff;
	// Switches
	int nHarmonics;
	int candidate_algorithm
	
	MSD_Parameters(){
		sigma_cutoff = 6.0;
		nHarmonics   = 8;
		candidate_algorithm = 0;
	}
};

#endif
