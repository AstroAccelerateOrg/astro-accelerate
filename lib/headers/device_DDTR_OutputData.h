#ifndef __ASTROACCELERATE_DDTR_OUTPUTDATA__
#define __ASTROACCELERATE_DDTR_OUTPUTDATA__

class DDTR_InputData {
public:
	int nRanges;
	int *inBin;
	int *nDMs_per_range;
	int *nSamples_per_range;
	
	float ***output_dmt_data;
	
	size_t output_in_bytes;
	
	//unsigned short *input_data;
	
	void Init(){
		inputsize_in_bytes = nsamp*nchans*sizeof(unsigned short);
	}
	
	void Allocate_host_input(){
		input_data = (unsigned short *) malloc(inputsize_in_bytes);
	}
	
	DDTR_InputData(){
		nchans = 0;
		nsamples = 0;
		nsamp = 0;
		nifs = 0;
		nbits = 0;
		tsamp = 0;
		tstart = 0;
		fch1 = 0;
		foff = 0;
		inputsize_in_bytes = 0;
		//input_data = NULL;
	}
	
	~DDTR_InputData(){
		//if(input_data!=NULL) free(input_data);
	}
};

#endif
