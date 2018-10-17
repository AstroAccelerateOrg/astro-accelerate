#ifndef ASTRO_ACCELERATE_SPS_SEARCH_HPP
#define ASTRO_ACCELERATE_SPS_SEARCH_HPP

#include <tuple>

#include "device_SPS_plan.hpp"
#include "device_SPS_parameters.hpp"
#include "device_MSD_parameters.hpp"
#include "device_SPS_DataDescription.hpp"

#include "device_analysis.hpp"
// TODO: Treat this as a proper search plan
class SPS_Search {
private:
	float *h_input; // I'm not sure if I want to include input from the host
	float *d_input;
	
	// NOTE: Moved from SPS Parameters - I want this to be the property of the whole search
	bool verbose;

    SPS_Plan spsplan;

protected:

	size_t max_candidates;
	size_t number_candidates;
	float *h_candidate_list;

public:
	
	SPS_Search(){
		h_input = NULL;
		d_input = NULL;
		h_candidate_list = NULL;
	}

	~SPS_Search(){
		if(h_candidate_list!=NULL) free(h_candidate_list);
	}

	/**
	 * @brief Check the verbosity mode set for the single pulse search
	 * 
	 * @return true verbose mode is on
	 * @return false verbose mode is off
	 */
	// TODO: think about including different verbosity levels - do we need all the information all the time?
	bool CheckVerbosity(void) {
		return verbose;
	}

	/**
	 * @brief Clears the candidate list
	 * 
	 */
	void ClearCandidateList(void){
		if (h_candidate_list != NULL) {
			free(h_candidate_list);
			h_candidate_list = NULL;
		}
		number_candidates = 0;
	}

    /**
     * @brief Creates a single pulse search plan
     * 
     */
    void CreateSearchPlan() {
        spsplan.Setup();
    }

	/**
	 * @brief Wrapper around the PrintSPSPlan method
	 * 
	 */
	// TODO: There must be a better way of doing that
	void PrintSearchPlan(void) const {
		spsplan.PrintSPSPlan();
	}

	/**
	 * @brief Starts the single pulse search
	 * 
	 * @param d_input pointer to the device memory containing dedispersed time series
	 * @return int 
	 */
	int RunSearch(float *d_input) {
		if (d_input == NULL) {
			std::cerr << "ERROR: input data pointer is NULL!" << std::endl;
			return(1);
		}
		number_candidates = 0;
		max_candidates = spsplan.GetMaxCandidates();
		h_candidate_list = (float*) malloc(max_candidates * 4 * sizeof(float));
		if (h_candidate_list == NULL) {
			std::cerr << "ERROR: not enough memory to allocate candidate list" << std::endl;
			return(1);
		}			
		
		analysis_GPU(verbose, d_input, h_candidate_list, number_candidates, spsplan);
		
		return(0);
	}

	/**
	 * @brief Saves the single pulse search candidates to file
	 * 
	 *  We still need to agree on the output data format and how to propagate it between different pipeline stages
	 * 
	 * @return int No error checking at the moment anyway
	 */
    int export_SPSData(void) const {
		char filename[200];
		
		/*
		for(int f=0; f<number_candidates; f++){
			h_candidate_list[4*list_pos]   = h_candidate_list[4*list_pos]*dm_step + dm_low;
			h_candidate_list[4*list_pos+1] = h_candidate_list[4*list_pos+1]*sampling_time + start_time;
			h_candidate_list[4*list_pos+2] = h_candidate_list[4*list_pos+2];
			h_candidate_list[4*list_pos+3] = h_candidate_list[4*list_pos+3]*inBin;
		}
		*/
		
		if (spsplan.GetSPSAlgorithm() == 0){
			sprintf(filename, "peak_analysed-t_%.2f-dm_%.2f-%.2f.dat", spsplan.GetStartTime(), std::get<0>(spsplan.GetDMLimits()), std::get<1>(spsplan.GetDMLimits()));
		}
		else if (spsplan.GetSPSAlgorithm() == 1) {
			sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", spsplan.GetStartTime(), std::get<0>(spsplan.GetDMLimits()), std::get<1>(spsplan.GetDMLimits()));
		}
					
		FILE *fp_out;
		
		if(number_candidates>0){
			if (( fp_out = fopen(filename, "wb") ) == NULL) return(1);
			fwrite(h_candidate_list, number_candidates*sizeof(float), 4, fp_out);
			fclose(fp_out);
		}
		return(0);
	}
	
	/* CandidateSubListPointer exportToSubList(void){
		CandidateSubListPointer newclp = new SPS_CandidateSubList(number_candidates, 0, h_candidate_list, NULL, NULL);
		newclp->time_start    = SPS_data.time_start;
		newclp->sampling_time = SPS_data.sampling_time;
		newclp->dm_step       = SPS_data.dm_step;
		newclp->dm_low        = SPS_data.dm_low;
		newclp->dm_high       = SPS_data.dm_high;
		newclp->inBin         = SPS_data.inBin;
		return(newclp);
	}
	*/
};

#endif