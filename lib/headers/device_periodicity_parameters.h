#ifndef __PER_PARAM__
#define __PER_PARAM__

class Periodicity_parameters {
public:
	float sigma_cutoff;
	int nHarmonics;
	int candidate_algorithm;
	int enable_outlier_rejection;
	float bln_sigma_constant;
	int export_powers;
	
	void assign(float t_sigma_cutoff, int t_nHarmonics, int t_candidate_algorithm, int t_enable_outlier_rejection, float t_bln_sigma_constant, int t_export_powers){
		sigma_cutoff             = t_sigma_cutoff;
		nHarmonics               = t_nHarmonics;
		candidate_algorithm      = t_candidate_algorithm;
		enable_outlier_rejection = t_enable_outlier_rejection;
		bln_sigma_constant       = t_bln_sigma_constant;
		export_powers            = t_export_powers;
	}
	
	void print_parameters(void){
		printf("Periodicity search parameters:\n");
		printf("Sigma cutoff: %f\n", sigma_cutoff);
		printf("Number of harmonics: %d\n", nHarmonics);
		printf("Candidate selection: %d\n", candidate_algorithm);
		printf("MSD variant: %s\n", (enable_outlier_rejection==1?"Baseline":"Normal") );
		printf("MSD outlier elimination sigma: %f\n", bln_sigma_constant);
		printf("Export powers:       %d\n", export_powers);
	}
};

#endif
