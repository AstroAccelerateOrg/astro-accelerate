#ifndef __ASTROACCELERATE_SPS_DATA_DESCRIPTION__
#define __ASTROACCELERATE_SPS_DATA_DESCRIPTION__

class SPS_Data_Description { // t-DM data
public:
	// Coordinate tranformation indexes-> [time, DM]
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
	
	SPS_Data_Description(float t_time_start, float t_sampling_time, float t_dm_step, float t_dm_low, float t_dm_high, int t_inBin, size_t t_nTimesamples, size_t t_nDMs){
		time_start = t_time_start;
		sampling_time = t_sampling_time;
		dm_step = t_dm_step;
		dm_low  = t_dm_low;
		dm_high = t_dm_high;
		inBin = t_inBin;
		nTimesamples = t_nTimesamples;
		nDMs = t_nDMs;
	}
};

#endif