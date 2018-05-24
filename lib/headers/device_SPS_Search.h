#ifndef __ASTROACCELERATE_SPS_SEARCH__
#define __ASTROACCELERATE_SPS_SEARCH__



class SPS_Search {
private:
	SPS_DataDescription SPS_data;
	SPS_Parameters      SPS_params;

protected:
	float *h_peaks;

public:

SPS_Search(){
	h_peaks = NULL;
}

~SPS_Search(){
	if(h_peaks!=NULL) free(h_peaks);
}

};

#endif
