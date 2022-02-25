#ifndef ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP
#define ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP

#include <stdio.h>
#include "aa_dedispersion_range.hpp"
#include "aa_ddtr_strategy.hpp"

#define PSR_HRMS_GREEDY 1
#define PSR_HRMS_PRESTO_PLUS 2
#define PSR_HRMS_PRESTO 3

namespace astroaccelerate {

/**
 * \class aa_periodicity_plan aa_periodicity_plan.hpp "include/aa_periodicity_plan.hpp"
 * \brief Class to set a periodicity plan.
 * \details A periodicity plan is required to create a periodicity strategy.
 * \author AstroAccelerate team.
 */
class aa_periodicity_plan {
public:
	/** \brief Trivial constructor for aa_periodicity_plan. */
	aa_periodicity_plan() {
		c_sigma_cutoff = 10.0;
		c_enable_outlier_rejection = true;
		c_sigma_outlier_rejection_threshold = 3.0;
		
		c_nHarmonics = 32;
		c_harmonic_sum_algorithm = 1;
		
		c_enable_interbinning = true;
		c_enable_scalloping_loss_mitigation = true;
		c_enable_spectrum_whitening = true;
		c_pad_to_nearest_higher_pow2 = false;
		
		c_candidate_selection_algorithm = 0;
	}
	
	/** \brief Constructor for aa_periodicity_plan that sets important member data on construction. */
	aa_periodicity_plan(
		const aa_ddtr_strategy &ddtr_strategy,
		float t_sigma_cutoff,
		bool t_enable_outlier_rejection,
		float t_sigma_outlier_rejection_threshold,
		int t_nHarmonics,
		int t_candidate_selection_algorithm
	) {
		c_sigma_cutoff = t_sigma_cutoff;
		c_enable_outlier_rejection = t_enable_outlier_rejection;
		c_sigma_outlier_rejection_threshold = t_sigma_outlier_rejection_threshold;
		c_nHarmonics = t_nHarmonics;
		c_harmonic_sum_algorithm = 1;
		c_candidate_selection_algorithm = t_candidate_selection_algorithm;
		
		c_enable_scalloping_loss_mitigation = true;
		c_enable_interbinning = true;
		c_enable_spectrum_whitening = true;
		c_pad_to_nearest_higher_pow2 = false;
		
		Extract_data_from_ddtr_strategy(ddtr_strategy);
	}
	
	/** \brief Constructor for aa_periodicity_plan that sets all member data on construction. */
	aa_periodicity_plan(
		const aa_ddtr_strategy &ddtr_strategy,
		float t_sigma_cutoff,
		bool t_enable_outlier_rejection,
		float t_sigma_outlier_rejection_threshold,
		int t_nHarmonics,
		int t_harmonic_sum_algorithm,
		int t_candidate_selection_algorithm,
		bool t_enable_scalloping_loss_mitigation,
		bool t_enable_interbinning,
		bool t_enable_spectrum_whitening
	) {
		c_sigma_cutoff = t_sigma_cutoff;
		c_enable_outlier_rejection = t_enable_outlier_rejection;
		c_sigma_outlier_rejection_threshold = t_sigma_outlier_rejection_threshold;
		c_nHarmonics = t_nHarmonics;
		c_harmonic_sum_algorithm = t_harmonic_sum_algorithm;
		c_candidate_selection_algorithm = t_candidate_selection_algorithm;
		
		c_enable_scalloping_loss_mitigation = t_enable_scalloping_loss_mitigation;
		c_enable_interbinning = t_enable_interbinning;
		c_enable_spectrum_whitening = t_enable_spectrum_whitening;
		c_pad_to_nearest_higher_pow2 = false;
		
		Extract_data_from_ddtr_strategy(ddtr_strategy);
	}
	
	/** \brief Constructor for aa_periodicity_plan that sets all member data on construction. */
	aa_periodicity_plan(
		float t_sigma_cutoff,
		bool t_enable_outlier_rejection,
		float t_sigma_outlier_rejection_threshold,
		int t_nHarmonics,
		int t_harmonic_sum_algorithm,
		
		bool t_enable_scalloping_loss_mitigation,
		bool t_enable_interbinning,
		
		int t_candidate_selection_algorithm,
		
		size_t t_processed_time_samples,
		double t_sampling_time,
		int    t_nRanges,
		float  *t_dm_low,
		float  *t_dm_high,
		float  *t_dm_step,
		int    *t_inBin,
		int const*const ndms
	) {
		c_sigma_cutoff = t_sigma_cutoff;
		c_enable_outlier_rejection = t_enable_outlier_rejection;
		c_sigma_outlier_rejection_threshold = t_sigma_outlier_rejection_threshold;
		c_nHarmonics = t_nHarmonics;
		c_harmonic_sum_algorithm = t_harmonic_sum_algorithm;
		c_candidate_selection_algorithm = t_candidate_selection_algorithm;
		
		c_enable_scalloping_loss_mitigation = t_enable_scalloping_loss_mitigation;
		c_enable_interbinning = t_enable_interbinning;
		c_enable_spectrum_whitening = true;
		c_pad_to_nearest_higher_pow2 = false;
		
		size_t c_nTimesamples = t_processed_time_samples;
		double c_sampling_time = t_sampling_time;
		printf("Creating plan from periodicity.\n");
		printf("Number of total number of time samples: %zu;\n", c_nTimesamples);
		printf("Sampling time: %e;\n", c_sampling_time);
		int nRanges = t_nRanges;
		for (int f = 0; f < nRanges; f++) {
			double dm_low = (double) t_dm_low[f];
			double dm_high = (double) t_dm_high[f];
			double dm_step = (double) t_dm_step[f];
			int inbin = t_inBin[f];
			int nDMs = ndms[f];
			aa_dedispersion_range newrange(dm_low, dm_high, dm_step, inbin, c_nTimesamples, nDMs, c_sampling_time);
			printf("dm_low=%f; dm_high=%f; dm_step=%f; inbin=%d; nTimesamples=%zu; nDMs=%d; sampling_time=%e;\n", dm_low, dm_high, dm_step, inbin, c_nTimesamples, nDMs, c_sampling_time); 
			ddtr_ranges.push_back(newrange);
		}
	}

	//----------- Setters ------------>
	void set_sigma_cutoff(float t_sigma_cutoff){
		c_sigma_cutoff = t_sigma_cutoff;
	}
	
	void set_sigma_outlier_rejection_threshold(float t_sigma_outlier_rejection_threshold){
		c_sigma_outlier_rejection_threshold = t_sigma_outlier_rejection_threshold;
	}
	
	void set_nHarmonics(int t_nHarmonics){
		c_nHarmonics = t_nHarmonics;
	}
	
	void set_outlier_rejection(bool t_switch){
		c_enable_outlier_rejection = t_switch;
	}
	
	void set_interbinning(bool t_switch){
		c_enable_interbinning = t_switch;
	}
	
	void set_scalloping_loss_mitigation(bool t_switch){
		c_enable_scalloping_loss_mitigation = t_switch;
	}
	
	void set_spectrum_whitening(bool t_switch){
		c_enable_spectrum_whitening = t_switch;
	}
	
	bool set_candidate_selection_algorithm(int algorithm_type){
		if(algorithm_type==0) {
			c_candidate_selection_algorithm = algorithm_type; // peak find
			return(true);
		}
		if(algorithm_type==1) {
			c_candidate_selection_algorithm = algorithm_type; // threshold
			return(true);
		}
		if(algorithm_type==2) {
			c_candidate_selection_algorithm = algorithm_type; // peak filtering
			return(true);
		}
		c_candidate_selection_algorithm = 1;
		return(false);
	}
	
	bool set_harmonic_sum_algorithm(int algorithm_type){
		if(algorithm_type==PSR_HRMS_GREEDY){
			c_harmonic_sum_algorithm = PSR_HRMS_GREEDY;
			return(true);
		}
		if(algorithm_type==PSR_HRMS_PRESTO_PLUS){
			c_harmonic_sum_algorithm = PSR_HRMS_PRESTO_PLUS;
			return(true);
		}
		if(algorithm_type==PSR_HRMS_PRESTO){
			c_harmonic_sum_algorithm = PSR_HRMS_PRESTO;
			return(true);
		}
		c_harmonic_sum_algorithm = PSR_HRMS_GREEDY;
		return(false);
	}
	
	//----------- Getters ------------>
	float sigma_cutoff() {
		return(c_sigma_cutoff);
	}
	
	float sigma_outlier_rejection_threshold(){
		return(c_sigma_outlier_rejection_threshold);
	}
	
	int nHarmonics(){
		return(c_nHarmonics);
	}
	
	int harmonic_sum_algorithm(){
		return(c_harmonic_sum_algorithm);
	}
	
	int candidate_selection_algorithm(){
		return(c_candidate_selection_algorithm);
	}
	
	bool enable_outlier_rejection(){
		return(c_enable_outlier_rejection);
	}
	
	bool enable_interbinning(){
		return(c_enable_interbinning);
	}
	
	bool enable_scalloping_loss_mitigation(){
		return(c_enable_scalloping_loss_mitigation);
	}
	
	bool enable_spectrum_whitening(){
		return(c_enable_spectrum_whitening);
	}
	
	bool pad_to_nearest_higher_pow2(){
		return(c_pad_to_nearest_higher_pow2);
	}
	
	size_t nRanges() const {
		return(ddtr_ranges.size());
	}
	
	aa_dedispersion_range get_range(size_t id) const {
		return(ddtr_ranges.at(id));
	}
	
	std::vector<aa_dedispersion_range> get_all_ranges(){
		return ddtr_ranges;
	}
	
	double sampling_time(){
		return(c_sampling_time);
	}
	
	void print(){
		printf("sigma_cutoff = %f;\n", c_sigma_cutoff);
		printf("sigma_outlier_rejection_threshold = %f;\n", c_sigma_outlier_rejection_threshold);
		printf("nHarmonics = %d;\n", c_nHarmonics);
		if(c_harmonic_sum_algorithm==0) printf("harmonic sum algorithm = simple;\n");
		else if(c_harmonic_sum_algorithm==1) printf("harmonic sum algorithm = greedy;\n");
		else if(c_harmonic_sum_algorithm==2) printf("harmonic sum algorithm = presto;\n");
		printf("candidate_selection_algorithm = %d;\n", c_candidate_selection_algorithm);
		printf("enable_outlier_rejection = "); if(c_enable_outlier_rejection) printf("true;\n"); else printf("false;\n");
		printf("enable_interbinning = "); if(c_enable_interbinning) printf("true;\n"); else printf("false;\n");
		printf("enable_scalloping_loss_mitigation = "); if(c_enable_scalloping_loss_mitigation) printf("true;\n"); else printf("false;\n");
		printf("enable_spectrum_whitening = "); if(c_enable_spectrum_whitening) printf("true;\n"); else printf("false;\n");
		printf("pad_to_nearest_higher_pow2 = "); if(c_pad_to_nearest_higher_pow2) printf("true;\n"); else printf("false;\n");
	}
	
	//------------------------- PRIVATE ------------------------------------
private:
	float c_sigma_cutoff;                      /** User selected value of sigma cutoff. Any event with SNR below this value will not be selected as candidate */
	float c_sigma_outlier_rejection_threshold; /** User selected sigma for outlier rejection. Any value with sigma greater than this will be rejected. */
	int   c_nHarmonics;                        /** User selected number of harmonics performed by the periodicity search. */
	int   c_harmonic_sum_algorithm;            /** Flag for selection of the harmonic sum algorithm */
	int   c_candidate_selection_algorithm;     /** User selected flag to select reduction algorithm to select candidates. */
	bool  c_enable_outlier_rejection;          /** User selected flag to enable or disable outlier rejection when calculating mean and standard deviation. */
	bool  c_enable_interbinning;               /** Enable or disable interpolation for periodicity search */
	bool  c_enable_scalloping_loss_mitigation; /** Enable 5-point convolution which mitigates scalloping loss */
	bool  c_enable_spectrum_whitening;         /** Enable or disable spectrum whitening (removal of the red noise) for periodicity search */
	bool  c_pad_to_nearest_higher_pow2;        /** Whether the periodicity will pad data to the nearest power of two for the Fourier transform. True by default */

	std::vector<aa_dedispersion_range> ddtr_ranges; /** DDTR - plan as calculated by ddtr strategy. */
	size_t c_nTimesamples;
	double c_sampling_time;
	
	void Extract_data_from_ddtr_strategy(const aa_ddtr_strategy &ddtr_strategy){
		//-------- Extraction from ddtr_strategy
		c_nTimesamples = 0;
		int nTimechunks = ddtr_strategy.num_tchunks();
		const std::vector<std::vector<int>> processed_samples = ddtr_strategy.t_processed();
		for(int f = 0; f < nTimechunks; f++){
			c_nTimesamples = c_nTimesamples + processed_samples[0][f];
		}
		printf("Extracting data from ddtr_strategy!\n");
		printf("Number of total number of time samples: %zu;\n", c_nTimesamples);
		
		aa_filterbank_metadata temp_metadata;
		temp_metadata = ddtr_strategy.metadata();
		double c_sampling_time = temp_metadata.tsamp();
		printf("Sampling time: %e;\n", c_sampling_time);
		
		int nRanges = ddtr_strategy.ndms_size();
		for (int f = 0; f < nRanges; f++) {
			aa_ddtr_plan::dm DM_range = ddtr_strategy.dm(f);
			double dm_low = (double) DM_range.low;
			double dm_high = (double) DM_range.high;
			double dm_step = (double) DM_range.step;
			int inbin = DM_range.inBin;
			int nDMs = ddtr_strategy.ndms(f);
			aa_dedispersion_range newrange(dm_low, dm_high, dm_step, inbin, c_nTimesamples, nDMs, c_sampling_time);
			printf("dm_low=%f; dm_high=%f; dm_step=%f; inbin=%d; nTimesamples=%zu; nDMs=%d; sampling_time=%e;\n", dm_low, dm_high, dm_step, inbin, c_nTimesamples, nDMs, c_sampling_time); 
			ddtr_ranges.push_back(newrange);
		}
	}

};
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP
