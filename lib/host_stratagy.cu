#include "headers/params.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "headers/headers_mains.h"

void stratagy(DDTR_Plan *DDTR_plan, size_t gpu_memory, DDTR_InputData *DDTR_data, AA_Parameters *AA_params){
	// This method relies on defining points when nsamps is a multiple of
	// nchans - bin on the diagonal or a fraction of it.
	
	double SPDT_fraction = 3.0/4.0; // 1.0 for MSD plane profile validation
	
	DDTR_plan->nchans = DDTR_data->nchans;
	DDTR_plan->nsamp  = DDTR_data->nsamp;
	DDTR_plan->nbits  = DDTR_data->nbits;
	DDTR_plan->tsamp  = DDTR_data->tsamp;
	
	
	// defining temporary variables
	float power    = DDTR_plan->power; // defined in constructor of DDTR_Plan
	float tsamp    = DDTR_plan->tsamp; 
	float fch1     = DDTR_data->fch1;
	float foff     = DDTR_data->foff;
	size_t nchans  = DDTR_data->nchans;
	size_t nsamp   = DDTR_data->nsamp;
	int nRanges    = DDTR_plan->nRanges;
	int i, j, c;
	int maxshift_high = 0;
	float n;
	float fmin = ( fch1 + ( foff*nchans ) );
	float fmin_pow = powf(fmin, power);
	float fmax_pow = powf(fch1, power);
	
	int error = 0;
	error = DDTR_plan->Allocate_ddtr_ranges(); if(error>0) exit(99);
	error = DDTR_plan->Allocate_dmshifts(); if(error>0) exit(99);
	//*dm_low  = (float *) malloc(( range ) * sizeof(float));
	//*dm_high = (float *) malloc(( range ) * sizeof(float));
	//*dm_step = (float *) malloc(( range ) * sizeof(float));
	//*ndms    = (int *) malloc(( range ) * sizeof(int));

	//*dmshifts = (float *) malloc(nchans * sizeof(float));
	

	//{{{ Calculate maxshift, the number of dms for this bin and
	//the highest value of dm to be calculated in this bin

	if (power != 2.0) {
		// Calculate time independent dm shifts
		for (c = 0; c < nchans; c++) {
			DDTR_plan->dmshifts[c] = 4148.741601f * ( ( 1.0 / pow(( fch1 + ( foff * c ) ), power) ) - ( 1.0 / pow(fch1, power) ) );
		}
	}
	else {
		// Calculate time independent dm shifts
		for (c = 0; c < nchans; c++) {
			DDTR_plan->dmshifts[c] = (float) ( 4148.741601f * ( ( 1.0 / pow((double) ( fch1 + ( foff * c ) ), power) ) - ( 1.0 / pow((double) fch1, power) ) ) );
		}
	}

	for (i = 0; i < nRanges; i++)	{
		modff(      ( ( (int) ( ( DDTR_plan->user_dm_high[i] - DDTR_plan->user_dm_low[i] ) / DDTR_plan->user_dm_step[i] ) + SDIVINDM ) / SDIVINDM )     , &n); // This calculates number of SDIVINDM blocks per DM range
		DDTR_plan->ndms[i] = (int) ( (int) n*SDIVINDM ); // This is number of DM trial per DM range
		if (DDTR_plan->max_ndms < DDTR_plan->ndms[i])
			DDTR_plan->max_ndms = DDTR_plan->ndms[i]; // looking for maximum number of DM trials for memory allocation
		DDTR_plan->total_ndms = DDTR_plan->total_ndms + DDTR_plan->ndms[i];
	}
	if(AA_params->verbose) printf("\nMaximum number of dm trials in any of the range steps:\t%d", DDTR_plan->max_ndms);

	DDTR_plan->dm_low[0] = DDTR_plan->user_dm_low[0];                        // 
	DDTR_plan->dm_high[0] = DDTR_plan->dm_low[0] + ( DDTR_plan->ndms[0]*DDTR_plan->user_dm_step[0] );   // Redefines DM plan to suit GPU
	DDTR_plan->dm_step[0] = DDTR_plan->user_dm_step[0];                      // 
	for (i = 1; i < nRanges; i++)	{
		DDTR_plan->dm_low[i] = DDTR_plan->dm_high[i - 1];
		DDTR_plan->dm_high[i] = DDTR_plan->dm_low[i] + DDTR_plan->ndms[i]*DDTR_plan->user_dm_step[i];
		DDTR_plan->dm_step[i] = DDTR_plan->user_dm_step[i];

		if (DDTR_plan->inBin[i - 1] > 1) {
			DDTR_plan->maxshift = (int) ceil(( ( DDTR_plan->dm_low[i - 1] + DDTR_plan->dm_step[i - 1]*DDTR_plan->ndms[i - 1] )*DDTR_plan->dmshifts[nchans - 1] ) / ( tsamp ));
			DDTR_plan->maxshift = (int) ceil((float) ( DDTR_plan->maxshift + ( (float) ( SDIVINT*2*SNUMREG ) ) ) / (float) DDTR_plan->inBin[i - 1]) / (float) ( SDIVINT*2*SNUMREG );
			DDTR_plan->maxshift = DDTR_plan->maxshift*( SDIVINT*2*SNUMREG )*DDTR_plan->inBin[i - 1];
			if ( DDTR_plan->maxshift > maxshift_high)
				maxshift_high = DDTR_plan->maxshift;
		}
	}

	if (DDTR_plan->inBin[nRanges - 1] > 1) {
		DDTR_plan->maxshift = (int) ceil(( ( DDTR_plan->dm_low[nRanges - 1] + DDTR_plan->dm_step[nRanges - 1] * DDTR_plan->ndms[nRanges - 1] ) * DDTR_plan->dmshifts[nchans - 1] ) / ( tsamp ));
		DDTR_plan->maxshift = (int) ceil((float) ( DDTR_plan->maxshift + ( (float) ( SDIVINT*2*SNUMREG ) ) ) / (float) DDTR_plan->inBin[nRanges - 1]) / (float) ( SDIVINT*2*SNUMREG );
		DDTR_plan->maxshift = DDTR_plan->maxshift * ( SDIVINT*2*SNUMREG )*DDTR_plan->inBin[nRanges - 1];
		if( DDTR_plan->maxshift > maxshift_high)
			maxshift_high = DDTR_plan->maxshift;
	}

	if (maxshift_high == 0)	{
		maxshift_high = (int) ceil(( ( DDTR_plan->dm_low[nRanges - 1] + DDTR_plan->dm_step[nRanges - 1] * ( DDTR_plan->ndms[nRanges - 1] ) ) * DDTR_plan->dmshifts[nchans - 1] ) / tsamp);
	}
	DDTR_plan->max_dm = ceil( DDTR_plan->dm_high[nRanges - 1] );

	DDTR_plan->maxshift = ( maxshift_high + ( SNUMREG * 2 * SDIVINT ) );
	printf("\nRange:\t%d, MAXSHIFT:\t%d, Scrunch value:\t%d", nRanges - 1, DDTR_plan->maxshift, DDTR_plan->inBin[nRanges - 1]);
	printf("\nMaximum dispersive delay:\t%.2f (s)", DDTR_plan->maxshift*tsamp);

	if (DDTR_plan->maxshift >= nsamp)	{
		printf("\n\nERROR!! Your maximum DM trial exceeds the number of samples you have.\nReduce your maximum DM trial\n\n");
		exit(1);
	}

	printf("\nDiagonal DM:\t%f", ( tsamp * nchans * 0.0001205 * powf(( fch1 + ( foff * ( nchans / 2 ) ) ), 3.0) ) / ( -foff * nchans ));
	if (DDTR_plan->maxshift >= nsamp)	{
		printf("ERROR!! Your maximum DM trial exceeds the number of samples you have.\nReduce your maximum DM trial");
		exit(1);
	}

	/* Four cases:
	 * 1) nchans < max_ndms & nsamp fits in GPU RAM
	 * 2) nchans > max_ndms & nsamp fits in GPU RAM
	 * 3) nchans < max_ndms & nsamp does not fit in GPU RAM
	 * 4) nchans > max_ndms & nsamp does not fit in GPU RAM
	 */

	unsigned int max_tsamps;

	// Allocate memory to store the t_processed ranges:
	DDTR_plan->Allocate_t_processed_outer();
	//t_processed = (int **) malloc(range * sizeof(int *));

	if (nchans < DDTR_plan->max_ndms) {
		// This means that we can cornerturn into the allocated output buffer 
		// without increasing the memory needed

		// Maximum number of samples we can fit in our GPU RAM is then given by:
		//max_tsamps = (unsigned int) ( (*gpu_memory) / ( sizeof(unsigned short) * ( (*max_ndms) + nchans ) ) ); // maximum number of timesamples we can fit into GPU memory
		size_t SPDT_memory_requirements = (AA_params->enable_analysis==1 ? (sizeof(float)*DDTR_plan->max_ndms*SPDT_fraction) : 0 );
		max_tsamps = (unsigned int) ( gpu_memory / ( sizeof(unsigned short)*nchans + sizeof(float)*DDTR_plan->max_ndms + SPDT_memory_requirements )); // maximum number of timesamples we can fit into GPU memory
		
		// Check that we dont have an out of range maxshift:
		if ( DDTR_plan->maxshift > max_tsamps)	{
			printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
			printf("\nReduce your maximum dm or increase the size of your dm step");
			exit(0);
		}

		// Next check to see if nsamp fits in GPU RAM:
		if (nsamp < max_tsamps)	{
			// We have case 1)
			// Allocate memory to hold the values of nsamps to be processed
			unsigned long int local_t_processed = (unsigned long int) floor(( (double) ( nsamp - DDTR_plan->maxshift ) / (double) DDTR_plan->inBin[nRanges - 1] ) / (double) ( SDIVINT*2*SNUMREG )); //number of timesamples per block
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * DDTR_plan->inBin[nRanges - 1];
			DDTR_plan->Allocate_t_processed_inner(1);
			for (i = 0; i < nRanges; i++) {
				DDTR_plan->t_processed[i] = (int *) malloc(sizeof(int)); // TODO: change to size_t
				DDTR_plan->t_processed[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) DDTR_plan->inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				DDTR_plan->t_processed[i][0] = DDTR_plan->t_processed[i][0] * ( SDIVINT*2*SNUMREG );
			}
			//( *num_tchunks ) = 1;
			printf("\nIn 1\n");
		}
		else {
			// We have case 3)
			// Work out how many time samples we can fit into ram 
			int samp_block_size = max_tsamps - DDTR_plan->maxshift;

			// Work out how many blocks of time samples we need to complete the processing
			// upto nsamp-maxshift
			//int num_blocks = (int) floor(( (float) nsamp - ( *maxshift ) )) / ( (float) ( samp_block_size ) ) + 1;

			// Find the common integer amount of samples between all bins
			int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) DDTR_plan->inBin[nRanges - 1] ) / (float) ( SDIVINT*2*SNUMREG ));
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * DDTR_plan->inBin[nRanges - 1];
			
			int num_blocks = (int) floor(( (float) nsamp - (float)DDTR_plan->maxshift )) / ( (float) ( local_t_processed ) );

			// Work out the remaining fraction to be processed
			int remainder =  nsamp -  (num_blocks*local_t_processed ) - DDTR_plan->maxshift;
			remainder = (int) floor((float) remainder / (float) DDTR_plan->inBin[nRanges - 1]) / (float) ( SDIVINT*2*SNUMREG );
			remainder = remainder * ( SDIVINT*2*SNUMREG ) * DDTR_plan->inBin[nRanges - 1];

			DDTR_plan->Allocate_t_processed_inner(num_blocks+1);
			for (i = 0; i < nRanges; i++)	{
				// Allocate memory to hold the values of nsamps to be processed
				//( *t_processed )[i] = (int *) malloc((num_blocks + 1) * sizeof(int));
				// Remember the last block holds less!
				for (j = 0; j < num_blocks; j++) {
					DDTR_plan->t_processed[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) DDTR_plan->inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
					DDTR_plan->t_processed[i][j] = DDTR_plan->t_processed[i][j] * ( SDIVINT*2*SNUMREG );
				}
				// fractional bit
				DDTR_plan->t_processed[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) DDTR_plan->inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				DDTR_plan->t_processed[i][num_blocks] = DDTR_plan->t_processed[i][num_blocks] * ( SDIVINT*2*SNUMREG );
			}
			//( *num_tchunks ) = num_blocks + 1;
			printf("\nIn 3\n");
			printf("\nnum_blocks:\t%d", num_blocks);
		}
	}
	else {
		// This means that we cannot cornerturn into the allocated output buffer 
		// without increasing the memory needed. Set the output buffer to be as large as the input buffer:

		// Maximum number of samples we can fit in our GPU RAM is then given by:
		//max_tsamps = (unsigned int) ( ( *gpu_memory ) / ( nchans * ( sizeof(float) + 2 * sizeof(unsigned short) ) ) );
		size_t SPDT_memory_requirements = (AA_params->enable_analysis==1 ? (sizeof(float)*DDTR_plan->max_ndms*SPDT_fraction) : 0 );
		max_tsamps = (unsigned int) ( gpu_memory / ( nchans * ( sizeof(float) + sizeof(unsigned short) )+ SPDT_memory_requirements ));

		// Check that we dont have an out of range maxshift:
		if ( DDTR_plan->maxshift > max_tsamps) {
			printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
			printf("\nReduce your maximum dm or increase the size of your dm step");
			exit(0);
		}

		// Next check to see if nsamp fits in GPU RAM:
		if (nsamp < max_tsamps) {
			// We have case 2)
			// Allocate memory to hold the values of nsamps to be processed
			int local_t_processed = (int) floor(( (float) ( nsamp - DDTR_plan->maxshift ) / (float) DDTR_plan->inBin[nRanges - 1] ) / (float) ( SDIVINT*2*SNUMREG ));
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * DDTR_plan->inBin[nRanges - 1];
			DDTR_plan->Allocate_t_processed_inner(1);
			for (i = 0; i < nRanges; i++) {
				//( *t_processed )[i] = (int *) malloc(sizeof(int));
				DDTR_plan->t_processed[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) DDTR_plan->inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				DDTR_plan->t_processed[i][0] = DDTR_plan->t_processed[i][0] * ( SDIVINT*2*SNUMREG );
			}
			//( *num_tchunks ) = 1;
			printf("\nIn 2\n");
		}
		else {
			// We have case 4)
			// Work out how many time samples we can fit into ram 
			int samp_block_size = max_tsamps - DDTR_plan->maxshift;

			// Work out how many blocks of time samples we need to complete the processing
			// upto nsamp-maxshift
			//int num_blocks = (int) floor(( (float) nsamp - (float) ( *maxshift ) ) / ( (float) samp_block_size ));

			// Find the common integer amount of samples between all bins
			int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) DDTR_plan->inBin[nRanges - 1] ) / (float) ( SDIVINT*2*SNUMREG ));
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * DDTR_plan->inBin[nRanges - 1];
			
			// samp_block_size was not used to calculate remainder instead there is local_t_processed which might be different
			int num_blocks = (int) floor(( (float) nsamp - (float) DDTR_plan->maxshift ) / ( (float) local_t_processed ));

			// Work out the remaining fraction to be processed
			int remainder = nsamp - ( num_blocks * local_t_processed ) - DDTR_plan->maxshift;
			remainder = (int) floor((float) remainder / (float) DDTR_plan->inBin[nRanges - 1]) / (float) ( SDIVINT*2*SNUMREG );
			remainder = remainder * ( SDIVINT*2*SNUMREG ) * DDTR_plan->inBin[nRanges - 1];

			DDTR_plan->Allocate_t_processed_inner(num_blocks+1);
			for (i = 0; i < nRanges; i++)	{
				// Allocate memory to hold the values of nsamps to be processed
				//( *t_processed )[i] = (int *) malloc(( num_blocks + 1 ) * sizeof(int));
				// Remember the last block holds less!
				for (j = 0; j < num_blocks; j++) {
					DDTR_plan->t_processed[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) DDTR_plan->inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
					DDTR_plan->t_processed[i][j] = DDTR_plan->t_processed[i][j] * ( SDIVINT*2*SNUMREG );
				}
				// fractional bit
				DDTR_plan->t_processed[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) DDTR_plan->inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				DDTR_plan->t_processed[i][num_blocks] = DDTR_plan->t_processed[i][num_blocks] * ( SDIVINT*2*SNUMREG );
			}
			//( *num_tchunks ) = num_blocks + 1;
			printf("\nIn 4\n");
		}
	}
	printf("\nMaxshift memory needed:\t%lu MB", nchans * DDTR_plan->maxshift * sizeof(unsigned short) / 1024 / 1024);
	if (nchans < DDTR_plan->max_ndms ) {
		printf("\nOutput memory needed:\t%lu MB", DDTR_plan->max_ndms * DDTR_plan->maxshift * sizeof(float) / 1024 / 1024);
	}
	else {
		printf("\nOutput memory needed:\t%lu MB", nchans * DDTR_plan->maxshift * sizeof(float) / 1024 / 1024);
	}

}
