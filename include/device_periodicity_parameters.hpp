#ifndef ASTRO_ACCELERATE_DEVICE_PERIODICITY_PARAMETERS_HPP
#define ASTRO_ACCELERATE_DEVICE_PERIODICITY_PARAMETERS_HPP

namespace astroaccelerate {

class Periodicity_parameters {
public:
	float sigma_cutoff;
	int nHarmonics;
	int candidate_algorithm;
	int enable_outlier_rejection;
	float OR_sigma_multiplier;
	int export_powers;
	
	void assign(float t_sigma_cutoff, int t_nHarmonics, int t_candidate_algorithm, int t_enable_outlier_rejection, float t_OR_sigma_multiplier, int t_export_powers){
		sigma_cutoff             = t_sigma_cutoff;
		nHarmonics               = t_nHarmonics;
		candidate_algorithm      = t_candidate_algorithm;
		enable_outlier_rejection = t_enable_outlier_rejection;
		OR_sigma_multiplier      = t_OR_sigma_multiplier;
		export_powers            = t_export_powers;
	}
	
	void print_parameters(void){
		printf("Periodicity search parameters:\n");
		printf("Sigma cutoff: %f\n", sigma_cutoff);
		printf("Number of harmonics: %d\n", nHarmonics);
		printf("Candidate selection: %d\n", candidate_algorithm);
		printf("MSD variant: %s\n", (enable_outlier_rejection==1?"Baseline":"Normal") );
		printf("MSD outlier elimination sigma: %f\n", OR_sigma_multiplier);
		printf("Export powers:       %d\n", export_powers);
	}
};

} //namespace astroaccelerate

#endif
