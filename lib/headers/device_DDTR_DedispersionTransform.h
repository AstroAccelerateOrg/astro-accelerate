#ifndef __ASTROACCELERATE_DDTR_DEDISPERSION__
#define __ASTROACCELERATE_DDTR_DEDISPERSION__

#include "headers/device_AA_Parameters.h"
#include "headers/device_DDTR_Plan.h"

//#include "headers/device_load_data.h"


class DDTR_DedispersionTransform {
private:
	// Bookkeeping: stores where we are in de-dispersion plan
	int old_timechunk, old_range;// initialize to -1? 
	int new_timechunk, new_range;// initialize to -1?
	size_t processed_samples; //inc = number of processed time samples so far.

public:

	DDTR_DedispersionTransform(){
		old_range = -1;
		new_range = -1;
		old_timechunk = -1;
		new_timechunk = -1;
		processed_samples = 0;
	}
	
	~DDTR_DedispersionTransform(){
		
	}

	size_t getProcessedSamples(){
		return(processed_samples);
	}
	
	int getCurrentTimeChunk(){
		return(new_timechunk);
	}
	
	int getCurrentRange(){
		return(new_range);
	}
	
	
	int endOfDDTRPlan( DDTR_Plan *DDTR_plan ){
		int l_cur_timechunk, l_cur_range;
		
		// getting index of the next time chunk and range and if we should decimate
		if(old_timechunk==-1 && old_range==-1 && DDTR_plan->nRanges!=0 && DDTR_plan->num_tchunks!=0){
			printf("endOfDDTRPlan: beginning of the DDTR\n");
			return(0);
		}
		else {
			l_cur_range = new_range;
			l_cur_range++;
			printf("endOfDDTRPlan: range: %d ", l_cur_range);
			if(l_cur_range==DDTR_plan->nRanges){
				l_cur_timechunk = new_timechunk;
				l_cur_timechunk++;
				printf("and new time-chunk: %d", l_cur_timechunk);
			}
			printf("\n");
		}
		
		// Exit check when all data are processed
		if(l_cur_timechunk==DDTR_plan->num_tchunks) {
			printf("endOfDDTRPlan: end of DDTR\n");
			return(1);
		}
		printf("endOfDDTRPlan: range: %d; time_chunk: %d\n", l_cur_range, l_cur_timechunk);
		return(0);
	}


	int dedisperseNextChunk(float *d_DDTR_output, unsigned short *d_DDTR_input, unsigned short  *h_input_buffer, DDTR_Plan *DDTR_plan, AA_Parameters *AA_params){

		// getting index of the next time chunk and range and if we should decimate
		if(old_timechunk==-1 && old_range==-1){
			// Initialize everything this is first run
			printf("DDTR.dedisperseNextChunk: Beginning\n");
			old_timechunk = -1;
			new_timechunk = 0;
			old_range = 0;
			new_range = 0;
		}
		else {
			old_timechunk = new_timechunk;
			old_range = new_range;
			new_range++;
			if(new_range==DDTR_plan->nRanges){
				old_timechunk = new_timechunk;
				new_timechunk++;
				old_range = 0;
				new_range = 0;
				printf("DDTR.dedisperseNextChunk: New timechunk\n");
				processed_samples = processed_samples + DDTR_plan->t_processed[0][old_timechunk];
			}
		}
		printf("DDTR.dedisperseNextChunk: new_range: %d; new time chunk:%d || old range: %d; old time chunk: %d\n",new_range, new_timechunk, old_range, old_timechunk);
		
		// Exit check when all data are processed
		if(new_timechunk==DDTR_plan->num_tchunks) return(-1);
		
		// Check if next de-dispersion step is correct
		if( new_range==0 && DDTR_plan->inBin[new_range]>1 ) {
			printf("Error: first range is not in full time resolution (has inBin>1)!\n");
			return(1);
		}
		if( DDTR_plan->inBin[new_range] > 2*DDTR_plan->inBin[old_range] ) {
			printf("Error: time resolution decreased too much between ranges %d and %d!\n", old_range, new_range);
			return(1);
		}
		
		// Set local variables which depend on inBin
		size_t local_maxshift = DDTR_plan->maxshift/DDTR_plan->inBin[new_range];
		float  local_tsamp    = DDTR_plan->tsamp*DDTR_plan->inBin[new_range];
		
		// Set other local variables
		float  tstart_local   = DDTR_plan->tsamp*processed_samples;
		
		// Initiate new time-chunk
		if(new_timechunk != old_timechunk){
			printf("t_processed: %d; time chunk: %d; maxshift: %d;\n", DDTR_plan->t_processed[0][new_timechunk], new_timechunk, local_maxshift);
			
			// This should be split into two functions
			load_data(-1, DDTR_plan->inBin, d_DDTR_input, &h_input_buffer[processed_samples*DDTR_plan->nchans], (int) DDTR_plan->t_processed[0][new_timechunk], local_maxshift, (int) DDTR_plan->nchans, DDTR_plan->dmshifts);
			
			checkCudaErrors(cudaGetLastError());
			
			if (AA_params->enable_zero_dm) {
				// ideally do as a callback
				zero_dm(d_DDTR_input, (int) DDTR_plan->nchans, (int) (DDTR_plan->t_processed[0][new_timechunk]+local_maxshift));
			}
			
			checkCudaErrors(cudaGetLastError());
			
			if (AA_params->enable_zero_dm_with_outliers) {
				zero_dm_outliers(d_DDTR_input, (int) DDTR_plan->nchans, (int) (DDTR_plan->t_processed[0][new_timechunk]+local_maxshift));
			}
			
			checkCudaErrors(cudaGetLastError());
		
			corner_turn(d_DDTR_input, d_DDTR_output, (int) DDTR_plan->nchans, (int) (DDTR_plan->t_processed[0][new_timechunk] + local_maxshift));
			
			checkCudaErrors(cudaGetLastError());
			
			if (AA_params->enable_old_rfi) {
				printf("\nPerforming old GPU rfi...");
				rfi_gpu(d_DDTR_input, (int) DDTR_plan->nchans, (int)(DDTR_plan->t_processed[0][new_timechunk]+local_maxshift));
			}
			
			checkCudaErrors(cudaGetLastError());
		}
		
		// -----------------------> Start processing data
		// ----> Set up stuff
		printf("\n\n%f\t%f\t%f\t%d\n", DDTR_plan->dm_low[new_range], DDTR_plan->dm_high[new_range], DDTR_plan->dm_step[new_range], DDTR_plan->ndms[new_range]);
		printf("Amount of telescope time processed: %f\n", tstart_local);
		
		cudaDeviceSynchronize();
		
		load_data(new_range, DDTR_plan->inBin, d_DDTR_input, &h_input_buffer[processed_samples*DDTR_plan->nchans], (int) DDTR_plan->t_processed[new_range][new_timechunk], (int) local_maxshift, (int) DDTR_plan->nchans, DDTR_plan->dmshifts);
		
		checkCudaErrors(cudaGetLastError());
		
		if( DDTR_plan->inBin[new_range] == 2*DDTR_plan->inBin[old_range] ) {
			bin_gpu(d_DDTR_input, d_DDTR_output, (int) DDTR_plan->nchans, (int) (DDTR_plan->t_processed[new_range - 1][new_timechunk] + local_maxshift*DDTR_plan->inBin[new_range]));
		}
		
		checkCudaErrors(cudaGetLastError());
				
		dedisperse(new_range, (int) DDTR_plan->t_processed[new_range][new_timechunk], DDTR_plan->inBin, DDTR_plan->dmshifts, d_DDTR_input, d_DDTR_output, (int) DDTR_plan->nchans, (int) ( DDTR_plan->t_processed[new_range][new_timechunk] + local_maxshift ), (int) local_maxshift, &local_tsamp, DDTR_plan->dm_low, DDTR_plan->dm_high, DDTR_plan->dm_step, DDTR_plan->ndms, DDTR_plan->nbits, AA_params->failsafe);
			
		checkCudaErrors(cudaGetLastError());
		
		return(0);
	}
};

#endif