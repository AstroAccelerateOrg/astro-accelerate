#ifndef ASTRO_ACCELERATE_SPS_PLAN_HPP
#define ASTRO_ACCELERATE_SPS_PLAN_HPP

#include <iostream>
#include <tuple>
#include <vector>

#include "device_MSD_parameters.hpp"
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
        int max_iterations;
        int output_shift;
        int shift;
        int start_taps;
        int unprocessed_samples;
        // NOTE: Replaced with total_unprocessed
        //int total_ut;
        int total_unprocessed;

        size_t max_candidates;
        size_t max_boxcar_width_desired;
        size_t max_boxcar_width_performed;
        // NOTE: Pulled from SPS Parameters
	    // Constants
	    float max_boxcar_width_in_sec;
	    float sigma_cutoff;
	    size_t max_candidate_array_size_in_bytes;
        
	    // Switches
	    int candidate_algorithm;
        std::vector<int> decimate_bc_limits;
        std::vector<int> boxcar_widths;

        // TODO: Keep MSD as a separate module - new options are likely to be added, so safet to keep it independent
        // TODO: Just store a local copy of the MSD module for SPS Plan
        // NOTE: Pulled from MSD Parameters
        float OR_sigma_multiplier;
	    // Switches
	    int enable_outlier_rejection;

        FilterBankDataMeta dataprop;
        MSD_Parameters msdparams;

    protected:

    public:

    public:
        SPS_Plan(void)
            : candidate_algorithm(0)
            , decimated_timesamples(0)
            , dtm(0)
            , iteration(0)
            , max_boxcar_width_in_sec(0.0f)
            , max_boxcar_width_desired(0)
            , max_boxcar_width_performed(0)
            , max_candidates(0)
            , max_candidate_array_size_in_bytes(0)
            , number_boxcars(0)
            , number_blocks(0)
            , output_shift(0)
            , shift(0)
            , sigma_cutoff(0.0f)
            , start_taps(0)
            , total_unprocessed(0)
            , unprocessed_samples(0) {

            dataprop.dm_low = 0.0f;
            dataprop.dm_high = 0.0f;
            dataprop.dm_step = 0.0f;
            dataprop.number_dms = 0;

            dataprop.binning_factor = 0;
            dataprop.start_time = 0.0f;
            dataprop.sampling_time = 0.0f;
            dataprop.timesamples = 0;

            dataprop.binned_sampling_time = dataprop.sampling_time * dataprop.binning_factor; // = 0.0f at the startup


        }

        ~SPS_Plan(void) = default;

        void Setup() {
            max_candidates = static_cast<size_t>((dataprop. number_dms) * (dataprop.timesamples) * 0.25);
            max_boxcar_width_desired = static_cast<int>(max_boxcar_width_in_sec / (dataprop.sampling_time));

            // NOTE: Set maximum boxcar widths performed between decimations
            decimate_bc_limits.push_back(PD_MAXTAPS);   
            decimate_bc_limits.push_back(16);
            decimate_bc_limits.push_back(16);
            decimate_bc_limits.push_back(16);
            decimate_bc_limits.push_back(8);  
        }

        /**
         * @brief Prints our basic single pulse search parameters
         */
        void PrintSPSPlan(void) const {
            std::cout << "Current Single Pulse Search Plan:" << std::endl;
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

        size_t GetCurrentMaxBoxcarWidth(void) const {
            return max_boxcar_width_performed;
        }

        std::tuple<float, float, float> GetDMLimits(void) const {
            return std::make_tuple(dataprop.dm_low, dataprop.dm_high, dataprop.dm_step);
        }

        /**
         * @brief Returns boxcar width for a given widths vector and element within the vector.
         * 
         * The function checks for the out-of-range access - returns the last element if element beyond the size of the vector is requested.
         * 
         * @param element index of the element to return 
         * @param bc_vector vector of boxcar widths to read values from
         * @return int the requested boxcar width
         */
    	int GetBCWidth(int element) const {
            // NOTE: There is enough elements in the widths list
            if (element < BC_widths.size()) {
                return BC_widths.at(element);
            // NOTE: Not enough elements - get the last element
            } else if (BC_widths.size() > 0) {
                return BC_widths.back();
            // NOTE: Edge case - the widths vector is empty
            } else {
                return 0;
            }
	    }

        size_t GetMaxCandidates(void) const {
            return max_candidates;
        }

        size_t GetNumberDMs(void) const {
            return dataprop.number_dms;
        }

        size_t GetCurrentTimeSamples(void) const {
            return dataprop.timesamples;
        }

        int GetCurrentBinningFactor(void) const {
            return dataprop.binning_factor;
        }

        float GetCurrentSamplingTime(void) const {
            return dataprop.binned_sampling_time;
        }

        float GetCurrentStartTime(void) const {
            return dataprop.start_time;
        }

        float GetOriginalSamplingTime(void) const {
            return dataprop.sampling_time;
        }
        /**
         * @brief Get the maximum iteration
         * 
         * @param max_width_performed 
         * @param max_boxcar_width 
         * @return the number of averaging iterations that have to be run in orded to reach the desired boxcar width
         */
        int CalculateMaxIteration(int max_boxcar_width) {
            int iteration = 0;
            while (max_boxcar_width_performed < max_boxcar_width) {
                max_boxcar_width_performed += GetBCWidth(iteration) *  (1 << iteration);
                iteration++;
            }
            max_iterations = iteration;
            return iteration;
        }

        int GetMaxIteration(void) const {
            return max_iterations;
        }

        void CreateSPSplan(void){
            int elements_per_block;
            int nRest;

            if(max_boxcar_width_desired > dataprop.timesamples) 
                max_iterations = CalculateMaxIteration(dataprop.timesamples);
            else
                max_iterations = CalculateMaxIteration(max_boxcar_width_desired);
            
            if(max_iterations > 0) {
                // NOTE: Set in the constructor
                // PDmp.shift        = 0;
                // PDmp.output_shift = 0;
                // PDmp.startTaps    = 0;
                // PDmp.iteration    = 0;
                
                // PDmp.decimated_timesamples = nTimesamples;
                // PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
                // PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
                
                decimated_timesamples = dataprop.timesamples;
                // NOTE: Each iteration decimates the data by a factor of 2 in the time dimension
                dtm = (dataprop.timesamples) >> (iteration + 1);
                // NOTE: That creates an even number of decimated time samples
                dtm = dtm - (dtm & 1);

                // NOTE: That doesn't sound right really
                // PDmp.nBoxcars = get_BC_width(0);
                // Elements_per_block = PD_NTHREADS*2 - PDmp.nBoxcars;
                // itemp = PDmp.decimated_timesamples;
                // PDmp.nBlocks = itemp/Elements_per_block;
                // nRest = itemp - PDmp.nBlocks*Elements_per_block;
                // if(nRest>0) PDmp.nBlocks++;
                // PDmp.unprocessed_samples = PDmp.nBoxcars + 6; // 6 is from where?
                // if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
                // PDmp.total_ut = PDmp.unprocessed_samples;
                
                // TODO: What's the logic behind this?
                number_boxcars = GetBCWidth(0);
                elements_per_block = PD_NTHREADS * 2 - number_boxcars;
                number_blocks = decimated_timesamples / elements_per_block;
                nRest = decimated_timesamples - number_blocks * elements_per_block;

                if (nRest > 0) {
                    number_blocks++;
                }

                // TODO: What's the logic behind this equation?
                unprocessed_samples = number_boxcars + 6;

                if (decimated_timesamples < unprocessed_samples) {
                    number_blocks = 0;
                }

                total_unprocessed = unprocessed_samples;
            }
        }
        
        /**
         * @brief Updates the plan depending on the decimation iteration
         * 
         * @param iteration 
         */
        void UpdateSPSPlan(int iteration) {
            int elements_per_block;
            int nRest;

            // for(int f=1; f<max_iterations; f++){
            //     // These are based on previous values of PDmp
            //     // PDmp.shift        = PDmp.nBoxcars/2;
            //     // PDmp.output_shift = PDmp.output_shift + PDmp.decimated_timesamples;
            //     // PDmp.startTaps    = PDmp.startTaps + PDmp.nBoxcars*(1<<PDmp.iteration);
            //     // PDmp.iteration    = PDmp.iteration + 1;
                
            //     // Definition of new PDmp values
            //     PDmp.decimated_timesamples = PDmp.dtm;
            //     PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
            //     PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
                
            //     PDmp.nBoxcars = get_BC_width(f);



            //     Elements_per_block=PD_NTHREADS*2 - PDmp.nBoxcars;
            //     itemp = PDmp.decimated_timesamples;
            //     PDmp.nBlocks = itemp/Elements_per_block;
            //     nRest = itemp - PDmp.nBlocks*Elements_per_block;
            //     if(nRest>0) PDmp.nBlocks++;
            //     PDmp.unprocessed_samples = PDmp.unprocessed_samples/2 + PDmp.nBoxcars + 6;
            //     if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
            //     PDmp.total_ut = PDmp.unprocessed_samples*(1<<PDmp.iteration);
                
            //     PD_plan->push_back(PDmp);
            // }

            shift = number_boxcars / 2;
            output_shift = output_shift + decimated_timesamples;
            start_taps = start_taps + number_boxcars * (1 << iteration);

            decimated_timesamples = dtm;
            dtm = (dataprop.timesamples) >> (iteration + 1);
            dtm = dtm - (dtm & 1);

            // TODO: This has to access the original vector containing 32, 16, 16, 16, 8 rather than the ones generated properly
            number_boxcars = GetBCWidth(iteration);
            elements_per_block = PD_NTHREADS * 2 - number_boxcars;
            number_blocks = decimated_timesamples / elements_per_block;
            nRest = decimated_timesamples - number_blocks * elements_per_block;
            if (nRest > 0) {
                number_blocks++;
            }

            unprocessed_samples = unprocessed_samples / 2 + number_boxcars + 6;
            if (decimated_timesamples < unprocessed_samples) {
                number_blocks = 0;
            }

            total_unprocessed = unprocessed_samples * (1 << iteration);

        }

        std::vector<int> CreateListOfBoxcars(int max_width_performed){
            int current_decimation = 1;
            int width = 0;
            int f = 0;

            while (width < max_width_performed) {
                for (int b = 0; b < GetBCWidth(f); b++) {
                    width = width + current_decimation;
                    boxcar_widths.push_back(width);
                }
                current_decimation *= 2;
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
            return dataprop.start_time;
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