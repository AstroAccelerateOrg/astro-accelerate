#ifndef __ASTROACCELERATE_SPS_CANDIDATELIST__
#define __ASTROACCELERATE_SPS_CANDIDATELIST__

struct Candidate {
	float DM;
	float time; //!!
	float SNR;
	float width;
};

class SPS_CandidateList { // time-DM data type
private:
	float *candidates;
	float *MSD;
	
public:
	// Coordinate transformation [indexes] -> [time, DM]
	float time_start;
	float sampling_time;
	float dm_step;
	float dm_low;
	float dm_high;
	// Amount of time decimation
	int inBin;
	
	void with_time_start(float value){
		time_start = value;
	}
	
	
	SPS_DataDescription(float time_start, float t_sampling_time, float t_dm_step, float t_dm_low, float t_dm_high, int t_inBin, size_t t_nTimesamples, size_t t_nDMs){
		this.time_start    = time_start;
		sampling_time = t_sampling_time;
		dm_step       = t_dm_step;
		dm_low        = t_dm_low;
		dm_high       = t_dm_high;
		inBin         = t_inBin;
	}
	
	SPS_DataDescription(){
		time_start    = 0;
		sampling_time = 0;
		dm_step       = 0;
		dm_low        = 0;
		dm_high       = 0;
		inBin         = 0;
	}
};

#endif
