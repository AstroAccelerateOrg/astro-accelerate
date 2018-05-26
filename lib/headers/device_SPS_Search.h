#ifndef __ASTROACCELERATE_SPS_SEARCH__
#define __ASTROACCELERATE_SPS_SEARCH__

#include "headers/device_SPS_DataDescription.h"
#include "headers/device_SPS_Parameters.h"
#include "headers/device_MSD_Parameters.h"
#include "headers/device_analysis.h"


class SPS_Search {
private:
	float *h_input; // I'm not sure if I want to include input from the host
	float *d_input;
	
protected:
	size_t max_nCandidates;
	size_t nCandidates;
	float *h_candidate_list;

public:
	SPS_DataDescription SPS_data;
	SPS_Parameters SPS_params;
	MSD_Parameters MSD_params;

	void setMSDParameters(MSD_Parameters t_MSD_params){
		MSD_params = t_MSD_params;
	}
	
	void setMSDParameters(MSD_Parameters *t_MSD_params){
		MSD_params = *t_MSD_params;
	}

	void setDataDescription(SPS_DataDescription t_SPS_data){
		SPS_data = t_SPS_data;
	}

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

	void setParameters(SPS_Parameters t_SPS_params){
		SPS_params = t_SPS_params;
	}
	
	void setParameters(SPS_Parameters *t_SPS_params){
		SPS_params = *t_SPS_params;
	}

	void setInputData(float *t_input){
		d_input = t_input;
	}

	void getCandidates(){
		
	}

	int search(){
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
	
	void clear(){
		if(h_candidate_list!=NULL) {
			free(h_candidate_list);
			h_candidate_list = NULL;
		}
		nCandidates = 0;
	}

	SPS_Search(){
		h_input = NULL;
		d_input = NULL;
		h_candidate_list = NULL;
	}

	~SPS_Search(){
		if(h_candidate_list!=NULL) free(h_candidate_list);
	}

};

#endif
