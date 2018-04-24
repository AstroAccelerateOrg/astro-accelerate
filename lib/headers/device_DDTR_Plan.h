#ifndef __ASTROACCELERATE_DDTR_PLAN__
#define __ASTROACCELERATE_DDTR_PLAN__

class DDTR_Plan {
public:
	int nRanges;
	float *user_dm_low;
	float *user_dm_high;
	float *user_dm_step;
	int *inBin;
	
	float *dm_low;
	float *dm_high;
	float *dm_step;
	int *ndms;
	
	float *dmshifts; //size of nchans
	int **t_processed;
	
	int maxshift;
	int num_tchunks;
	int max_ndms;
	int total_ndms;
	float max_dm;
	
	// Unknowns stuff
	float power;
	
	// Sizes
	size_t gpu_inputsize;
	size_t gpu_outputsize;
	size_t host_inputsize;
	size_t host_outputsize;
	
	
	// Description of input data
	int nchans; //number of frequency channels
	int nsamp; // number of timesamples
	float tsamp; //sampling time 
	
	int Allocate_user_ranges(){
		int error=0;
		user_dm_low   = (float *) malloc( nRanges*sizeof(float) );
		if(user_dm_low==NULL) error++;
		user_dm_high  = (float *) malloc( nRanges*sizeof(float) );
		if(user_dm_high==NULL) error++;
		user_dm_step  = (float *) malloc( nRanges*sizeof(float) );
		if(user_dm_step==NULL) error++;
		inBin         = (float *) malloc( nRanges*sizeof(float) );
		if(inBin==NULL) error++;
		return(error);
	}
	
	int Allocate_ddtr_ranges(){
		int error=0;
		dm_low   = (float *) malloc( nRanges*sizeof(float) );
		if(dm_low==NULL) error++;
		dm_high  = (float *) malloc( nRanges*sizeof(float) );
		if(dm_high==NULL) error++;
		dm_step  = (float *) malloc( nRanges*sizeof(float) );
		if(dm_step==NULL) error++;
		ndms     = (float *) malloc( nRanges*sizeof(float) );
		if(ndms==NULL) error++;
		return(error);
	}
	
	int Allocate_dmshifts(){
		dmshifts = (float *) malloc( nchans*sizeof(float) );
		if(dmshifts==NULL) return(1);
	}
	
	int Allocate_t_processed_outer(){
		t_processed = (int **) malloc( nRanges*sizeof(int *) );
		if(t_processed==NULL) return(1);		
	}
	
	int Allocate_t_processed_inner(int t_num_tchunks){
		int error = 0;
		num_tchunks = t_num_tchunks;
		for(int r=0; r<nRanges; r++){
			t_processed[r] = (int *) malloc( num_tchunks*sizeof(int) );
			if(t_processed[r] == NULL ) error++;
		}
		return(error);
	}
	
	
	DDTR_Plan(){
		nRanges = 0;
		
		user_dm_low  = NULL;
		user_dm_high = NULL;
		user_dm_step = NULL;
		inBin        = NULL;
		
		dm_low  = NULL;
		dm_high = NULL;
		dm_step = NULL;
		ndms    = NULL;
		
		dmshifts = NULL;
		
		num_tchunks = 0;
		t_processed = NULL;
		
		maxshift   = 0;
		max_ndms   = 0;
		total_ndms = 0;
		max_dm     = 0;
		
		power = 2.0;
		
		gpu_inputsize   = 0;
		gpu_outputsize  = 0;
		host_inputsize  = 0;
		host_outputsize = 0;
		
		nchans = 0;
		nsamp  = 0;
		tsamp  = 0;
	}
	
	~DDTR_Plan(){
		if(user_dm_low!=NULL) free(user_dm_low);
		if(user_dm_high!=NULL) free(user_dm_high);
		if(user_dm_step!=NULL) free(user_dm_step);
		if(inBin!=NULL) free(inBin);
		
		if(dm_low!=NULL) free(dm_low);
		if(dm_high!=NULL) free(dm_high);
		if(dm_step!=NULL) free(dm_step);
		if(ndms!=NULL) free(ndms);
		
		if(dmshifts!=NULL) free(dmshifts);
		
		for(int r=0; r<nRanges, r++){
			if(t_processed[r]!=NULL) free(t_processed[r]); 
		}
		if(t_processed!=NULL) free(t_processed);
	}
};

#endif
