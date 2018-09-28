#ifndef ASTRO_ACCELERATE_SPS_PLAN_HPP
#define ASTRO_ACCELERATE_SPS_PLAN_HPP

#include <tuple>

#include "FilterBankDataMeta.hpp"
#include "params.hpp"
/**
 * @brief Structure that holds the basic properties of the SPS input data
 * 
 */

// NOTE: That has been moved from device_SPS_DataDescription.hpp
struct DataProperties {
    float dm_low;
    float dm_high;
    float dm_step;
    float start_time;
    float sampling_time;

    float binned_sampling_time;

    int binning_factor;

    size_t number_dms;
    size_t timesamples;
};

/**
 * @brief Class which holds the SPS plan
 * 
 */
class SPS_Plan {

    private:
        int decimated_timesamples;
        int dtm;    // decimated time samples
        int iteration;
        int number_boxcars;
        int number_blocks;
        int output_shift;
        int shift;
        int start_taps;
        int unprocessed_samples;
        int total_ut;

        size_t max_candidates;
        size_t max_boxcar_width;
        // NOTE: Pulled from SPS Parameters
	    // Constants
	    float max_boxcar_width_in_sec;
	    float sigma_cutoff;
	    size_t max_candidate_array_size_in_bytes;
        
	    // Switches
	    int candidate_algorithm;
        std::vector<int> BC_widths;

        // TODO: Keep MSD as a separate module - new options are likely to be added, so safet to keep it independent
        // TODO: Just store a local copy of the MSD module for SPS Plan
        // NOTE: Pulled from MSD Parameters
        float OR_sigma_multiplier;
	    // Switches
	    int enable_outlier_rejection;

        FilterBankDataMeta *dataprop;
        MSD_Parameters msdparams;

    protected:

    public:

    public:
        SPS_Plan(void)
            : candidate_algorithm(0)
            , decimated_timesamples(0)
            // NOTE: What is dtm?
            // TODO: Come up with a better name
            , dtm(0)
            // NOTE: What is iteration?
            // TODO: Is there a better name for it?
            , iteration(0)
            , max_boxcar_width_in_sec(0.0f)
            , max_candidates(0)
            , max_candidate_array_size_in_bytes(0)
            , number_boxcars(0)
            , number_blocks(0)
            // TODO: What is the difference between output_shift and shift?
            , output_shift(0)
            , shift(0)
            , sigma_cutoff(0.0f)
            // TODO: What is TAPS?
            , start_taps(0)
            // NOTE: That relates to the fact that final time samples are not processed properly
            // TODO: What does this do?
            , total_ut(0) 
            , unprocessed_samples(0) {

            dataprop -> dm_low = 0.0f;
            dataprop -> dm_high = 0.0f;
            dataprop -> dm_step = 0.0f;
            dataprop -> number_dms = 0;

            dataprop -> binning_factor = 0;
            dataprop -> start_time = 0.0f;
            dataprop -> sampling_time = 0.0f;
            dataprop -> timesamples = 0;

            dataprop -> binned_sampling_time = dataprop -> sampling_time * dataprop -> binning_factor;

            // NOTE: These were set just before the SPS run and reset after every iteration
            // TODO: Move this to a setup method - constructor can't throw
            BC_widths.push_back(PD_MAXTAPS);   
            BC_widths.push_back(16);
            BC_widths.push_back(16);
            BC_widths.push_back(16);
            BC_widths.push_back(8);  

        }

        ~SPS_Plan(void) = default;

        void Setup() {
            max_candidates = static_cast<size_t>((dataprop -> number_dms) * (dataprop -> timesamples) * 0.25);
            max_boxcar_width = static_cast<int>(max_boxcar_width_in_sec / (dataprop -> sampling_time));
        }

        /**
         * @brief Prints our basic single pulse search parameters
         */
        void PrintSPSPlan(void) const {
            std::cout << "SPS Plan:" << std::endl;
            std::cout << "\tmax_boxcar_with_in_sec: " << max_boxcar_width_in_sec << std::endl;
            std::cout << "\tsigma_cutoff: " << sigma_cutoff << std::endl;
            std::cout << "\tmax_candidate_array_size_in_bytes: " << max_candidate_array_size_in_bytes << std::endl;

            // NOTE: Technically there are 2 atm, but there might be more in the future?
            if (candidate_algorithm == 1) {
                std::cout << "\tcandidate algorithm: threshold" << std::endl;
            } else if (candidate_algorithm == 0) {
                std::cout << "\tcandidate algorithm: peak-find" << std::endl;
            }
            std::cout << "\tBC widths: ";
            if (BC_widths.size() == 0) {
                std::cout << "not set" << std::endl;
            } else {
                for (std::vector<int>::const_iterator iwidth = BC_widths.begin(); iwidth != BC_widths.end(); ++iwidth) {
                    std::cout << *iwidth << " ";
                }
            }
            std::cout << std::endl;
            std::cout << "---------------------<" << std::endl;
        }

        /**
         * @brief Create a new plan based on the user input
         * 
         */
        void CreateSPSPlan() {

        }

        /**
         * @brief Update the plan based on the user input
         * 
         */
        // TODO: Check which parameters actually need updating
        void UpdateSPSPlan() {

        }

        size_t GetCurrentMaxBoxcarWidth(void) const {
            return max_boxcar_width;
        }

        std::tuple<float, float, float> GetDMLimits(void) const {
            return std::make_tuple(dataprop.dm_low, dataprop.dm_high, dataprop.dm_step);
        }

        size_t GetMaxCandidates(void) const {
            return max_candidates;
        }

        size_t GetNumberDMs(void) const {
            return dataprop -> number_dms;
        }

        size_t GetTimeSamples(void) const {
            return dataprop -> timesamples;
        }

        int GetCurrentBinningFactor(void) const {
            return dataprop -> binning_factor;
        }

        float GetCurrentSamplingTime(void) const {
            return dataprop -> binned_sampling_time;
        }

        float GetCurrentStartTime(void) const {
            return dataprop -> start_time;
        }

        float GetOriginalSamplingTime(void) const {
            return dataprop -> sampling_time;
        }
        /**
         * @brief Get the maximum iteration
         * 
         * @param max_width_performed 
         * @param max_boxcar_width 
         * @return int 
         */
        int GetMaxIteration(int *max_width_performed, int max_boxcar_width){
            int start_taps, iteration, f;
            
            start_taps = 0;
            iteration = 0;
            f=0;
            while(start_taps<max_boxcar_width){
                start_taps = start_taps + BC_widths[f]*(1<<f);
                f = f + 1;
            }
            
            iteration = f;
            *max_width_performed=start_taps;
            return(iteration);
        }

	void Create_PD_plan(std::vector<PulseDetection_plan> *PD_plan, int *max_width_performed, int max_boxcar_width, size_t const nTimesamples){
		int Elements_per_block, itemp, nRest;
		int max_iteration;
		PulseDetection_plan PDmp;
		
		if(max_boxcar_width>nTimesamples) 
			max_iteration = get_max_iteration(max_width_performed, nTimesamples);
		else
			max_iteration = get_max_iteration(max_width_performed, max_boxcar_width);
		
		if(max_iteration>0){
			PDmp.shift        = 0;
			PDmp.output_shift = 0;
			PDmp.startTaps    = 0;
			PDmp.iteration    = 0;
			
			PDmp.decimated_timesamples = nTimesamples;
			PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
			PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
			
			PDmp.nBoxcars = get_BC_width(0);
			Elements_per_block = PD_NTHREADS*2 - PDmp.nBoxcars;
			itemp = PDmp.decimated_timesamples;
			PDmp.nBlocks = itemp/Elements_per_block;
			nRest = itemp - PDmp.nBlocks*Elements_per_block;
			if(nRest>0) PDmp.nBlocks++;
			PDmp.unprocessed_samples = PDmp.nBoxcars + 6; // 6 is from where?
			if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
			PDmp.total_ut = PDmp.unprocessed_samples;
			
			
			PD_plan->push_back(PDmp);
			
			for(int f=1; f<max_iteration; f++){
				// These are based on previous values of PDmp
				PDmp.shift        = PDmp.nBoxcars/2;
				PDmp.output_shift = PDmp.output_shift + PDmp.decimated_timesamples;
				PDmp.startTaps    = PDmp.startTaps + PDmp.nBoxcars*(1<<PDmp.iteration);
				PDmp.iteration    = PDmp.iteration + 1;
				
				// Definition of new PDmp values
				PDmp.decimated_timesamples = PDmp.dtm;
				PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
				PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
				
				PDmp.nBoxcars = get_BC_width(f);
				Elements_per_block=PD_NTHREADS*2 - PDmp.nBoxcars;
				itemp = PDmp.decimated_timesamples;
				PDmp.nBlocks = itemp/Elements_per_block;
				nRest = itemp - PDmp.nBlocks*Elements_per_block;
				if(nRest>0) PDmp.nBlocks++;
				PDmp.unprocessed_samples = PDmp.unprocessed_samples/2 + PDmp.nBoxcars + 6;
				if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
				PDmp.total_ut = PDmp.unprocessed_samples*(1<<PDmp.iteration);
				
				PD_plan->push_back(PDmp);
			}
		}
	}
	
	void Create_list_of_boxcar_widths(std::vector<int> *boxcar_widths, int max_width_performed){
		int DIT_value, DIT_factor, width, f;
		DIT_value = 1;
		DIT_factor = 2;
		width = 0;
		f=0;
		while(width<max_width_performed){
			for(int b=0; b<get_BC_width(f); b++){
				width = width + DIT_value;
				boxcar_widths->push_back(width);
			}
			DIT_value = DIT_value*DIT_factor;
			f++;
		}
	}

        /**
         * @brief Set the algorithm number
         * @param algorithm user suppkied algorithm number: 0 - peak-find, 1 - threshold 
         */
        void SetAlgorithm(short algorithm) {
            candidate_algorithm = algorithm;
        }

        /**
         * @brief Get the sigma cutoff for the SPS algorithm
         * 
         * @return float 
         */
        float GetSigmaCutoff(void) const {
            return sigma_cutoff;
        }

        /**
         * @brief Return the single-pulse detection algorithm currently in use
         * @return short
         */
        short GetSPSAlgorithm(void) const {
            return candidate_algorithm;
        }

        float GetStartTime(void) const {
            return dataprop -> start_time;
        }

        /**
         * @brief Adds boxcar width
         * 
         * @param width width of boxcar function in time samples
         */
        void AddBCWidth(int width) {
            BC_widths.push_back(width);
        }

            /**
         * @brief Clears the boxcar width vector
         * 
         */
        void ClearBCWidths(void){
            BC_widths.clear();
        }
};

#endif