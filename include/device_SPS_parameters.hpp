#ifndef ASTRO_ACCELERATE_SPS_PARAMETERS_HPP
#define ASTRO_ACCELERATE_SPS_PARAMETERS_HPP

#include <iostream>
#include <vector>

#include "params.hpp"

class SPS_Parameters {

		// TODO: Move public variables to private and write proper setters and getters
	private:

    public:
		float max_boxcar_width_in_sec;
		float sigma_cutoff;
		// TODO: This variable name is a tad too long
		size_t max_candidate_array_size_in_bytes;

        // TODO: Should that be made a string? thres / peak is more understandable and memorable than 0 / 1 I think
		short candidate_algorithm;

		// NOTE: verbose is not used anywhere in here
		// TODO: check if it is used somewhere in the code
		int verbose;

		std::vector<int> BC_widths;
	public:

		SPS_Parameters(void)
			: max_boxcar_width_in_sec(0.5f)
			, sigma_cutoff(6.0f)
			, max_candidate_array_size_in_bytes(0)
			, candidate_algorithm(0)
			, verbose(0) {
		
            // NOTE: These were set just before the SPS run and reset after every iteration
            BC_widths.push_back(PD_MAXTAPS);   
            BC_widths.push_back(16);
            BC_widths.push_back(16);
            BC_widths.push_back(16);
            BC_widths.push_back(8);  

		}

		~SPS_Parameters(void) {

		}

		/**
		 * @brief Prints our basic single pulse search parameters
		 */
		void PrintSPSParameters(void) const {
			std::cout << "SPS Parameters:" << std::endl;
			std::cout << "\tmax_boxcar_with_in_sec: " << max_boxcar_width_in_sec << std::endl;
			std::cout << "\tsigma_cutoff: " << sigma_cutoff << std::endl;
			std::cout << "\tmax_candidate_array_size_in_bytes: " << max_candidate_array_size_in_bytes << std::endl;

			// NOTE: Technically there are 2 atm, but there might be more in the future?
			if (candidate_algorithm == 1) {
				std::cout << "\tcandidate algorithm: threshold" << std::endl;
			} else if (candidate_algorithm == 0) {
				std::cout << "\tcandidate algorithm: peak-find" << std::endl;
			}
			printf("\tverbose: %d;\n", verbose);
			std::cout << "\tverbose: " << verbose << std::endl;
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
	 * @brief Clears the boxcar with vector
	 * 
	 */
	void ClearBCWidths(void){
		BC_widths.clear();
	}

    /**
     * @brief Get the maximum iteration
     * 
     * @param max_width_performed 
     * @param max_boxcar_width 
     * @return int 
     */
	int GetMaxIteration(int *max_width_performed, int max_boxcar_width){
		int startTaps, iteration, f;
		
		startTaps = 0;
		iteration = 0;
		f=0;
		while(startTaps<max_boxcar_width){
			startTaps = startTaps + get_BC_width(f)*(1<<f);
			f = f + 1;
		}
		
		iteration = f;
		*max_width_performed=startTaps;
		return(iteration);
	}

	/**
	 * @brief Set the algorithm number
	 * @param algorithm user suppkied algorithm number: 0 - peak-find, 1 - threshold 
	 */
	void SetAlgorithm(short algorithm) {
		candidate_algorithm = algorithm;
	}

	/**
	 * @brief Return the single-pulse detection algorithm currently in use
	 * @return short
	 */
	short GetAlgorithms(void) const {
		return candidate_algorithm;
	}

	/**
	 * @brief Adds boxcar width
	 * 
	 * @param width width of boxcar function in time samples
	 */
    void AddBCWidth(int width) {
        BC_widths.push_back(width);
    }
};


#endif