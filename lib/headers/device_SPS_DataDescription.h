#ifndef __ASTROACCELERATE_SPS_DATADESCRIPTION__
#define __ASTROACCELERATE_SPS_DATADESCRIPTION__

class SPS_DataDescription { // time-DM data type
public:
	// Coordinate transformation [indexes] -> [time, DM]
	float time_start;
	float sampling_time;
	float dm_step;
	float dm_low;
	float dm_high;
	// Amount of time decimation
	int inBin;
	// Data dimensions
	size_t nTimesamples; //t_processed
	size_t nDMs;
	
	void set(float t_time_start, float t_sampling_time, float t_dm_step, float t_dm_low, float t_dm_high, int t_inBin, size_t t_nTimesamples, size_t t_nDMs){
		time_start    = t_time_start;
		sampling_time = t_sampling_time;
		dm_step       = t_dm_step;
		dm_low        = t_dm_low;
		dm_high       = t_dm_high;
		inBin         = t_inBin;
		nTimesamples  = t_nTimesamples;
		nDMs          = t_nDMs;
	}
	
	SPS_DataDescription(){
		time_start    = 0;
		sampling_time = 0;
		dm_step       = 0;
		dm_low        = 0;
		dm_high       = 0;
		inBin         = 0;
		nTimesamples  = 0;
		nDMs          = 0;
	}
};

#endif
