#ifndef __ASTROACCELERATE_DDTR_INPUTDATA__
#define __ASTROACCELERATE_DDTR_INPUTDATA__

//TODO: rename stuff:
// nchans -> nChannels or nFrequencyChannels
// nsamp -> nSamples or nTimesamples
// nbits -> nBits?
// tsamp -> sampling_time
// fch1 -> Should not be necessary ro rename. It should not be accessed from outside and it is in the fil header?
// foff -> ---||---
// TODO: Add array which will hold the data.
// TODO: Add methods for loading data from files

class DDTR_InputData {
public:
	size_t nchans;
	size_t nsamples; // this in number of samples from header. It is NOT used in the code. Why?
	size_t nsamp;
	size_t nifs;
	int nbits;
	float tsamp;
	float tstart;
	float fch1;
	float foff;
	
	size_t inputsize_in_bytes;
	
	//unsigned short *input_data;
	
	void Init(){
		inputsize_in_bytes = nsamp*nchans*sizeof(unsigned short);
	}
	
	//void Allocate_host_input(){
	//	input_data = (unsigned short *) malloc(inputsize_in_bytes);
	//}
	
	//void Input_from_fil_file(){
	//}
	
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
