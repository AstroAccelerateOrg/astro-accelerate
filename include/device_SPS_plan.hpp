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

struct ProcessingDetails {
        int decimated_timesamples;
        int dtm;
        int iteration;
        int number_blocks;
        int number_boxcars;
        int output_shift;
        int shift;
        int start_taps;
        int total_unprocessed;
        int unprocessed_samples;
};

class SPS_Plan {

    private:
        int max_iterations;
        // NOTE: Replaced with total_unprocessed
        //int total_ut;
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

        std::vector<ProcessingDetails> details;

    protected:

    public:

    public:
        SPS_Plan(void)
            : candidate_algorithm(0)
            , max_boxcar_width_in_sec(0.50f)
            , max_boxcar_width_desired(0)
            , max_boxcar_width_performed(0)
            , max_candidates(0)
            , max_candidate_array_size_in_bytes(0)
            , sigma_cutoff(6.0f) {

            // TODO: We still need to pass this information
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

        /**
         * @brief Returns boxcar width for a given widths vector and element within the vector.
         * 
         * The function checks for the out-of-range access - returns the last element if element beyond the size of the vector is requested.
         * 
         * @param element index of the element to return 
         * @param wodths vector of boxcar widths to read values from
         * @return int the requested boxcar width
         */
    	int GetBCWidth(int element, std::vector<int> &widths) const {
            // NOTE: There is enough elements in the widths list
            if (element < widths.size()) {
                return widths.at(element);
            // NOTE: Not enough elements - get the last element
            } else if (widths.size() > 0) {
                return widths.back();
            // NOTE: Edge case - the widths vector is empty
            } else {
                return 0;
            }
	    }
        
        /**
         * @brief Calculates the number of averaging iterations
         * 
         * @param max_width_performed 
         * @param max_boxcar_width 
         * @return the number of averaging iterations that have to be run in orded to reach the desired boxcar width
         */
        int CalculateMaxIteration(int max_boxcar_width) {
            int iteration = 0;
            max_boxcar_width_performed = 0;
            while (max_boxcar_width_performed < max_boxcar_width) {
                max_boxcar_width_performed += GetBCWidth(iteration, decimate_bc_limits) *  (1 << iteration);
                iteration++;
            }
            max_iterations = iteration;
            return iteration;
        }

        /**
         * @brief Creates a list of boxcar widths
         * 
         */
        void CreateListOfBoxcars(){
            int current_decimation = 1;
            int width = 0;
            int f = 0;

            while (width < max_boxcar_width_performed) {
                for (int b = 0; b < GetBCWidth(f, decimate_bc_limits); b++) {
                    width = width + current_decimation;
                    boxcar_widths.push_back(width);
                }
                current_decimation *= 2;
                f++;
            }
        }

        /**
         * @brief Creates the actual plan calculating all the necessary variables
         * 
         */
        // TODO: This should really accept a structure that contains all the necessary information about the current data chunk
        void Setup() {

            max_candidates = static_cast<size_t>((dataprop. number_dms) * (dataprop.timesamples) * 0.25);
            max_boxcar_width_desired = static_cast<int>(max_boxcar_width_in_sec / (dataprop.sampling_time));

            // NOTE: Set maximum boxcar widths performed between decimations
            decimate_bc_limits.push_back(PD_MAXTAPS);   
            decimate_bc_limits.push_back(16);
            decimate_bc_limits.push_back(16);
            decimate_bc_limits.push_back(16);
            decimate_bc_limits.push_back(8);  

            int elements_per_block;
            int elements_rem;

            ProcessingDetails tmpdetails = {0};

            if(max_boxcar_width_desired > dataprop.timesamples) 
                max_iterations = CalculateMaxIteration(dataprop.timesamples);
            else
                max_iterations = CalculateMaxIteration(max_boxcar_width_desired);
            
            if(max_iterations > 0) {
                tmpdetails.shift = 0;
                tmpdetails.output_shift = 0;
                tmpdetails.start_taps = 0;
                tmpdetails.iteration = 0;

                tmpdetails.decimated_timesamples = dataprop.timesamples;
                // NOTE: Each iteration decimates the data by a factor of 2 in the time dimension
                tmpdetails.dtm = (dataprop.timesamples) >> (tmpdetails.iteration + 1);
                // NOTE: That creates an even number of decimated time samples
                tmpdetails.dtm = tmpdetails.dtm - (tmpdetails.dtm & 1);
                
                // TODO: What's the logic behind this?
                tmpdetails.number_boxcars = GetBCWidth(0, decimate_bc_limits);
                elements_per_block = PD_NTHREADS * 2 - tmpdetails.number_boxcars;
                tmpdetails.number_blocks = tmpdetails.decimated_timesamples / elements_per_block;
                elements_rem = tmpdetails.decimated_timesamples - tmpdetails.number_blocks * elements_per_block;

                if (elements_rem > 0) {
                    tmpdetails.number_blocks++;
                }

                // TODO: What's the logic behind this equation?
                tmpdetails.unprocessed_samples = tmpdetails.number_boxcars + 6;

                if (tmpdetails.decimated_timesamples < tmpdetails.unprocessed_samples) {
                    tmpdetails.number_blocks = 0;
                }

                tmpdetails.total_unprocessed = tmpdetails.unprocessed_samples;

                details.push_back(tmpdetails);
            }

            for (int iiter = 0; iiter < max_iterations; ++iiter) {
                tmpdetails.shift = tmpdetails.number_boxcars / 2;
                tmpdetails.output_shift = tmpdetails.output_shift + tmpdetails.decimated_timesamples;
                tmpdetails.start_taps = tmpdetails.start_taps + tmpdetails.number_boxcars * (1 << tmpdetails.iteration);
                tmpdetails.iteration = tmpdetails.iteration + 1;

                tmpdetails.decimated_timesamples = tmpdetails.dtm;
                tmpdetails.dtm = (dataprop.timesamples) >> (tmpdetails.iteration + 1);
                tmpdetails.dtm = tmpdetails.dtm - (tmpdetails.dtm & 1);

                // TODO: This has to access the original vector containing 32, 16, 16, 16, 8 rather than the ones generated properly
                tmpdetails.number_boxcars = GetBCWidth(tmpdetails.iteration, decimate_bc_limits);
                elements_per_block = PD_NTHREADS * 2 - tmpdetails.number_boxcars;
                tmpdetails.number_blocks = tmpdetails.decimated_timesamples / elements_per_block;
                elements_rem = tmpdetails.decimated_timesamples - tmpdetails.number_blocks * elements_per_block;
                if (elements_rem > 0) {
                    tmpdetails.number_blocks++;
                }

                tmpdetails.unprocessed_samples = tmpdetails.unprocessed_samples / 2 + tmpdetails.number_boxcars + 6;
                if (tmpdetails.decimated_timesamples < tmpdetails.unprocessed_samples) {
                    tmpdetails.number_blocks = 0;
                }

                tmpdetails.total_unprocessed = tmpdetails.unprocessed_samples * (1 << tmpdetails.iteration);

                details.push_back(tmpdetails);
            }

            CreateListOfBoxcars();
        }

        /**
         * @brief Prints our basic single pulse search parameters
         */
        void PrintSPSPlan(void) const {
            std::cout << "Current Single Pulse Search Plan:" << std::endl;
            std::cout << "\tmax_boxcar_with_in_sec: " << max_boxcar_width_in_sec << std::endl;
            std::cout << "\tmax_boxcar_width_desired: " << max_boxcar_width_desired << std::endl;
            std::cout << "\tmax_boxcar_width_performed: " << max_boxcar_width_performed << std::endl;
            std::cout << "\tsigma_cutoff: " << sigma_cutoff << std::endl;
            std::cout << "\tmax_candidate_array_size_in_bytes: " << max_candidate_array_size_in_bytes << std::endl;

            // NOTE: Technically there are 2 atm, but there might be more in the future?
            if (candidate_algorithm == 1) {
                std::cout << "\tcandidate algorithm: threshold" << std::endl;
            } else if (candidate_algorithm == 0) {
                std::cout << "\tcandidate algorithm: peak-find" << std::endl;
            }

            std::cout << "\tBC widths: ";
            if (boxcar_widths.size() == 0) {
                std::cout << "not set" << std::endl;
            } else {
                for (std::vector<int>::const_iterator iwidth = boxcar_widths.begin(); iwidth != boxcar_widths.end(); ++iwidth) {
                    std::cout << *iwidth << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "\tNumber of decimation iterations: " << max_iterations << std::endl;
            if (decimate_bc_limits.size() > 0) {
                std::cout << "\tDecimation in time at the following BC widths: ";
                for (std::vector<int>::const_iterator iwidth = decimate_bc_limits.begin(); iwidth != decimate_bc_limits.end(); ++iwidth) {
                    std::cout << *iwidth << " ";
                }
            }

            std::cout << std::endl;
            std::cout << "---------------------<" << std::endl;
        }

        ProcessingDetails GetDetails(int iteration) {
            return details.at(iteration);
        }

        size_t GetCurrentMaxBoxcarWidth(void) const {
            return max_boxcar_width_performed;
        }

        std::tuple<float, float, float> GetDMLimits(void) const {
            return std::make_tuple(dataprop.dm_low, dataprop.dm_high, dataprop.dm_step);
        }
        /**
         * @brief Returns the number of candidates we can safely process
         * 
         * @return size_t 
         */
        size_t GetMaxCandidates(void) const {
            return max_candidates;
        }

        /**
         * @brief Returns the number of DM values in the current data chunk
         * 
         * @return size_t 
         */
        size_t GetNumberDMs(void) const {
            return dataprop.number_dms;
        }
        /**
         * @brief Returns the number of time samples in the current data chunk
         * 
         * @return size_t 
         */
        size_t GetCurrentTimeSamples(void) const {
            return dataprop.timesamples;
        }
        /**
         * @brief Returns the binning (averaging) factor for the current data chunk
         * 
         * @return int 
         */
        int GetCurrentBinningFactor(void) const {
            return dataprop.binning_factor;
        }

        /**
         * @brief Returns the current sampling time for the current data chunk
         * 
         * That includes the binning factors
         * 
         * @return float 
         */
        float GetCurrentSamplingTime(void) const {
            return dataprop.binned_sampling_time;
        }

        /**
         * @brief Returns the start time of the current data chunk
         * 
         * @return float 
         */
        float GetCurrentStartTime(void) const {
            return dataprop.start_time;
        }

        /**
         * @brief Returns the original sampling time of the current data chunk
         * 
         * This does not include any pre- and post-dedispersion averaging factors
         * 
         * @return float 
         */
        float GetOriginalSamplingTime(void) const {
            return dataprop.sampling_time;
        }

        /**
         * @brief Returns the number of time decimation iterations to be run during single pulse search
         * 
         * @return int 
         */
        int GetMaxIteration(void) const {
            return max_iterations;
        }

        std::vector<int> GetListOfBoxcars(void) const {
            return boxcar_widths;
        }

        /**
         * @brief Set the algorithm number
         * @param algorithm user supplied algorithm number: 0 - peak-find, 1 - threshold 
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

        MSD_Parameters GetMSDParameters(void) const {
            return msdparams;
        }

};

#endif