#ifndef __ASTROACCELERATE_DDTR_DATA__
#define __ASTROACCELERATE_DDTR_DATA__

class DDTR_Data {
public:
	size_t nchans;
	size_t nsamples;
	size_t nsamp;
	size_t nifs;
	int nbits;
	float tsamp;
	float tstart;
	float fch1;
	float foff;
	
	DDTR_Data(){
		nchans = 0;
		nsamples = 0;
		nsamp = 0;
		nifs = 0;
		nbits = 0;
		tsamp = 0;
		tstart = 0;
		fch1 = 0;
		foff = 0;
	}
};

#endif
