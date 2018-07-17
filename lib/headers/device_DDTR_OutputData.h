#ifndef __ASTROACCELERATE_DDTR_OUTPUTDATA__
#define __ASTROACCELERATE_DDTR_OUTPUTDATA__

// TODO: We could change dedispersion plan to something like in periodicity and input and output data could inherit properties from some kind of base class. Then it would be easier to go through dedispersion time chunks.
class DDTR_OutputData {
public:
	float  tstart_local;
	float  sampling_time;
	float  dm_step;
	float  dm_low;
	float  dm_high;
	int    inBin;
	size_t nTimesamples;
	int    nDMs;
	
	float *d_DDTR_output;
	
	DDTR_OutputData(){
		tstart_local = 0;
		sampling_time = 0;
		dm_step = 0;
		dm_low = 0;
		dm_high = 0;
		inBin = 0;
		nTimesamples = 0;
		nDMs = 0;
		d_DDTR_output = NULL;
	}
	
	~DDTR_OutputData(){
		// d_DDTR_output is only temporary pointer to the data. It should not be deallocated.
	}
};

#endif
