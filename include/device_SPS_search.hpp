#ifndef ASTRO_ACCELERATE_SPS_SEARCH_HPP
#define ASTRO_ACCELERATE_SPS_SEARCH_HPP

#include "device_SPS_plan.hpp"
#include "device_SPS_parameters.hpp"
#include "device_MSD_parameters.hpp"
#include "device_SPS_DataDescription.hpp"

// TODO: Treat this as a proper search plan
class SPS_Search {
private:
	float *h_input; // I'm not sure if I want to include input from the host
	float *d_input;
	
	// NOTE: Moved from SPS Parameters - I want this to be the property of the whole search
	bool verbose;

    SPS_Plan spsplan;

protected:

	size_t max_nCandidates;
	size_t nCandidates;
	float *h_candidate_list;

public:


	// NOTE: Now part of SPS_plan
	//SPS_DataDescription SPS_data;	
	//SPS_Parameters SPS_params;
	
	MSD_Parameters MSD_params;

	SPS_Search(){
		h_input = NULL;
		d_input = NULL;
		h_candidate_list = NULL;
	}

	~SPS_Search(){
		if(h_candidate_list!=NULL) free(h_candidate_list);
	}

	/**
	 * @brief 
	 * 
	 * @return true verbose mode is on
	 * @return false verbose mode is off
	 */
	// TODO: think about including different verbosity levels - do we need all the information all the time?
	bool CheckVerbosity(void) {
		return verbose;
	}

	void Clear(){
		if(h_candidate_list!=NULL) {
			free(h_candidate_list);
			h_candidate_list = NULL;
		}
		nCandidates = 0;
	}

	/**
	 * @brief Get a copy of the single-pulse search plan
	 * 
	 * @return SPS_Plan&
	 */
	SPS_Plan GetSPSPlan(void) {
		return spsplan;
	}

    /**
     * @brief Creates a single pulse search plan
     * 
     */
    void CreateSearchPlan() {
        spsplan.CreateSPSPlan()
    }

    /**
     * @brief Set the Data Description object
     * 
     * @param t_time_start 
     * @param t_sampling_time 
     * @param t_dm_step 
     * @param t_dm_low 
     * @param t_dm_high 
     * @param t_inBin 
     * @param t_nTimesamples 
     * @param t_nDMs 
     */
	void setDataDescription(float t_time_start, float t_sampling_time, float t_dm_step, float t_dm_low, float t_dm_high, int t_inBin, size_t t_nTimesamples, size_t t_nDMs){
		SPS_data.time_start = t_time_start;
		SPS_data.sampling_time = t_sampling_time;
		SPS_data.dm_step = t_dm_step;
		SPS_data.dm_low  = t_dm_low;
		SPS_data.dm_high = t_dm_high;
		SPS_data.inBin = t_inBin;
		SPS_data.nTimesamples = t_nTimesamples;
		SPS_data.nDMs = t_nDMs;
	}

	void SetSPSParameters(SPS_Parameters t_SPS_params){
		SPS_params = t_SPS_params;
	}

    /**
     * @brief Sets the pointer to the dedispersed time series
     * 
     * @param tmpinput pointer to the dedispersed time series memory
     */
	void SetSPSInputData(float *tmpinput){
		d_input = tmpinput;
	}

	void getCandidates() {
		
	}
	/**
	 * @brief Wrapper around the PrintSPSPlan method
	 * 
	 */
	// TODO: There must be a better way of doing that
	void PrintSearchPlan(void) {
		spsplan.PrintSPSPlan();
	}

	int RunSPS(void) {
		if(d_input==NULL) {
			if(SPS_params.verbose) printf("ERROR: input data pointer is NULL!\n");
			return(1);
		}
		max_nCandidates = 0;
		nCandidates = 0;
		max_nCandidates = (size_t) ( (SPS_data.nDMs*SPS_data.nTimesamples)/4 );
		h_candidate_list = (float*) malloc(max_nCandidates*4*sizeof(float));
		if(h_candidate_list==NULL) {
			if(SPS_params.verbose) printf("ERROR: not enough memory to allocate candidate list\n");
			return(1);
		}			
		
		analysis_GPU(h_candidate_list, &nCandidates, max_nCandidates, SPS_data, d_input, &SPS_params, &MSD_params);
		
		return(0);
	}
	
	void ClearCandidateList(){
		if(h_candidate_list!=NULL) {
			free(h_candidate_list);
			h_candidate_list = NULL;
		}
		nCandidates = 0;
	}

    int export_SPSData(void) const {
		char filename[200];
		
		/*
		for(int f=0; f<nCandidates; f++){
			h_candidate_list[4*list_pos]   = h_candidate_list[4*list_pos]*dm_step + dm_low;
			h_candidate_list[4*list_pos+1] = h_candidate_list[4*list_pos+1]*sampling_time + start_time;
			h_candidate_list[4*list_pos+2] = h_candidate_list[4*list_pos+2];
			h_candidate_list[4*list_pos+3] = h_candidate_list[4*list_pos+3]*inBin;
		}
		*/
		
		if(SPS_params.candidate_algorithm==1){
			sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", SPS_data.time_start, SPS_data.dm_low, SPS_data.dm_high);
		}
		else {
			sprintf(filename, "peak_analysed-t_%.2f-dm_%.2f-%.2f.dat", SPS_data.time_start, SPS_data.dm_low, SPS_data.dm_high);
		}
					
		FILE *fp_out;
		
		if(nCandidates>0){
			if (( fp_out = fopen(filename, "wb") ) == NULL) return(1);
			fwrite(h_candidate_list, nCandidates*sizeof(float), 4, fp_out);
			fclose(fp_out);
		}
		return(0);
	}
	
	CandidateSubListPointer exportToSubList(void){
		CandidateSubListPointer newclp = new SPS_CandidateSubList(nCandidates, 0, h_candidate_list, NULL, NULL);
		newclp->time_start    = SPS_data.time_start;
		newclp->sampling_time = SPS_data.sampling_time;
		newclp->dm_step       = SPS_data.dm_step;
		newclp->dm_low        = SPS_data.dm_low;
		newclp->dm_high       = SPS_data.dm_high;
		newclp->inBin         = SPS_data.inBin;
		return(newclp);
	}

};

#endif