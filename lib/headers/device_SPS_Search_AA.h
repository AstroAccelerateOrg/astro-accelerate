#ifndef __ASTROACCELERATE_SPS_SEARCH_AA__
#define __ASTROACCELERATE_SPS_SEARCH_AA__


#include "headers/device_SPS_Search.h"
#include "headers/device_SPS_DataDescription.h"
#include "headers/device_SPS_Parameters.h"
#include "headers/device_MSD_Parameters.h"
#include "headers/device_analysis.h"


class SPS_Search_AA: public SPS_Search {

public:
	int export_SPSData(){
		char filename[200];
		
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
};

#endif
