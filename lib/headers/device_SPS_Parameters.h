#ifndef __ASTROACCELERATE_SPS_PARAMETERS__
#define __ASTROACCELERATE_SPS_PARAMETERS__

#include <vector>
#include "headers/device_SPS_BC_plan.h"

class SPS_Parameters {
public:
	// Constants
	float max_boxcar_width_in_sec;
	float sigma_cutoff;
	
	// Switches
	int candidate_algorithm;
	int verbose;
	
	// SPS plan
	std::vector<int> BC_widths;
	
	SPS_Parameters(){
		max_boxcar_width_in_sec = 0.5f;
		sigma_cutoff = 6.0f;
		candidate_algorithm = 0;
		verbose = 0;
	}
	
	//----------------------> BC_widths
	void set_BC_widths(int *widths, int size){
		BC_widths.clear();
		for(int f=0; f<size; f++) BC_widths.push_back(widths[f]);
	}
	
	void set_BC_widths(std::vector<int> *widths){
		BC_widths.clear();
		for(int f=0; f<widths->size(); f++) BC_widths.push_back(widths->operator[](f));
	}
	
	void add_BC_width(int width){
		BC_widths.push_back(width);
	}
	
	int get_BC_width(int el){
		int last_element = BC_widths.size()-1;
		if( el>last_element ) {
			printf("Returning last element\n");
			return( BC_widths[last_element] );
		}
		else {
			printf("Returning normal element\n");
			return( BC_widths[el] );
		}
	}
	//-----------------------------------<
	
	int get_max_iteration(int *max_width_performed, int max_boxcar_width){
		int startTaps, iteration, f;
		
		startTaps = 0;
		iteration = 0;
		f=0;
		while(startTaps<max_boxcar_width){
			startTaps = startTaps + get_BC_width(f)*(1<<f);
			f = f + 1;
		}
		
		iteration = f;
		*max_width_performed=startTaps;
		return(iteration);
	}
	
	void Create_PD_plan(std::vector<PulseDetection_plan> *PD_plan, int *max_width_performed, int max_boxcar_width, size_t const nTimesamples){
		int Elements_per_block, itemp, nRest;
		int max_iteration;
		PulseDetection_plan PDmp;
		
		if(max_boxcar_width>nTimesamples) 
			max_iteration = get_max_iteration(max_width_performed, nTimesamples);
		else
			max_iteration = get_max_iteration(max_width_performed, max_boxcar_width);
		
		if(max_iteration>0){
			PDmp.shift        = 0;
			PDmp.output_shift = 0;
			PDmp.startTaps    = 0;
			PDmp.iteration    = 0;
			
			PDmp.decimated_timesamples = nTimesamples;
			PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
			PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
			
			PDmp.nBoxcars = get_BC_width(0);
			Elements_per_block = PD_NTHREADS*2 - PDmp.nBoxcars;
			itemp = PDmp.decimated_timesamples;
			PDmp.nBlocks = itemp/Elements_per_block;
			nRest = itemp - PDmp.nBlocks*Elements_per_block;
			if(nRest>0) PDmp.nBlocks++;
			PDmp.unprocessed_samples = PDmp.nBoxcars + 6; // 6 is from where?
			if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
			PDmp.total_ut = PDmp.unprocessed_samples;
			
			
			PD_plan->push_back(PDmp);
			
			for(int f=1; f<max_iteration; f++){
				// These are based on previous values of PDmp
				PDmp.shift        = PDmp.nBoxcars/2;
				PDmp.output_shift = PDmp.output_shift + PDmp.decimated_timesamples;
				PDmp.startTaps    = PDmp.startTaps + PDmp.nBoxcars*(1<<PDmp.iteration);
				PDmp.iteration    = PDmp.iteration + 1;
				
				// Definition of new PDmp values
				PDmp.decimated_timesamples = PDmp.dtm;
				PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
				PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
				
				PDmp.nBoxcars = get_BC_width(f);
				Elements_per_block=PD_NTHREADS*2 - PDmp.nBoxcars;
				itemp = PDmp.decimated_timesamples;
				PDmp.nBlocks = itemp/Elements_per_block;
				nRest = itemp - PDmp.nBlocks*Elements_per_block;
				if(nRest>0) PDmp.nBlocks++;
				PDmp.unprocessed_samples = PDmp.unprocessed_samples/2 + PDmp.nBoxcars + 6;
				if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
				PDmp.total_ut = PDmp.unprocessed_samples*(1<<PDmp.iteration);
				
				PD_plan->push_back(PDmp);
			}
		}
	}
};

#endif
