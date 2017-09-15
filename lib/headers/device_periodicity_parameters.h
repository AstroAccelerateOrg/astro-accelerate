#ifndef __PER_PARAM__
#define __PER_PARAM__

class Periodicity_parameters {
public:
	float sigma_cutoff;
	int nHarmonics;
	
	void assign(float t_sigma_cutoff, int t_nHarmonics){
		sigma_cutoff = t_sigma_cutoff;
		nHarmonics   = t_nHarmonics;
	}
	
	void print_parameters(void){
		printf("Periodicity search parameters:\n");
		printf("Sigma cutoff: %f\n", sigma_cutoff);
		printf("Number of harmonics: %d\n", nHarmonics);
	}
};

#endif
