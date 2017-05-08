#include "headers/params.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void stratagy(int *maxshift, int *max_samps, int *num_tchunks, int *max_ndms, int *total_ndms, float *max_dm, float power, int nchans, int nsamp, float fch1, float foff, float tsamp, int range, float *user_dm_low, float *user_dm_high, float *user_dm_step, float **dm_low, float **dm_high, float **dm_step, int **ndms, float **dmshifts, int *inBin, int ***t_processed, size_t *gpu_memory, size_t SPS_mem_requirement) {
	// This method relies on defining points when nsamps is a multiple of
	// nchans - bin on the diagonal or a fraction of it.

	int i, j, c;
	int maxshift_high = 0;

	float n;
	float fmin = ( fch1 + ( foff * nchans ) );
	float fmin_pow = powf(fmin, power);
	float fmax_pow = powf(fch1, power);

	*dm_low = (float *) malloc(( range ) * sizeof(float));
	*dm_high = (float *) malloc(( range ) * sizeof(float));
	*dm_step = (float *) malloc(( range ) * sizeof(float));
	*ndms = (int *) malloc(( range ) * sizeof(int));

	*dmshifts = (float *) malloc(nchans * sizeof(float));

	//{{{ Calculate maxshift, the number of dms for this bin and
	//the highest value of dm to be calculated in this bin

	if (power != 2.0) {
		// Calculate time independent dm shifts
		for (c = 0; c < nchans; c++) {
			( *dmshifts )[c] = 4148.741601f * ( ( 1.0 / pow(( fch1 + ( foff * c ) ), power) ) - ( 1.0 / pow(fch1, power) ) );
		}
	}
	else {
		// Calculate time independent dm shifts
		for (c = 0; c < nchans; c++) {
			( *dmshifts )[c] = (float) ( 4148.741601f * ( ( 1.0 / pow((double) ( fch1 + ( foff * c ) ), power) ) - ( 1.0 / pow((double) fch1, power) ) ) );
		}
	}

	for (i = 0; i < range; i++)	{
		modff(      ( ( (int) ( ( user_dm_high[i] - user_dm_low[i] ) / user_dm_step[i] ) + SDIVINDM ) / SDIVINDM )     , &n); // This calculates number of SDIVINDM blocks per DM range
		( *ndms )[i] = (int) ( (int) n * SDIVINDM ); // This is number of DM trial per DM range
		if (*max_ndms < ( *ndms )[i])
			*max_ndms = ( *ndms )[i]; // looking for maximum number of DM trials for memory allocation
		*total_ndms = *total_ndms + ( *ndms )[i];
	}
	printf("\nMaximum number of dm trials in any of the range steps:\t%d", *max_ndms);

	( *dm_low )[0] = user_dm_low[0];                        // 
	( *dm_high )[0] = ( *dm_low )[0] + ( ( *ndms )[0] * ( user_dm_step[0] ) );   // Redefines DM plan to suit GPU
	( *dm_step )[0] = user_dm_step[0];                      // 
	for (i = 1; i < range; i++)	{
		( *dm_low )[i] = ( *dm_high )[i - 1];
		( *dm_high )[i] = ( *dm_low )[i] + ( *ndms )[i] * user_dm_step[i];
		( *dm_step )[i] = user_dm_step[i];

		if (inBin[i - 1] > 1) {
			*maxshift = (int) ceil(( ( (*dm_low)[i - 1] + (*dm_step)[i - 1]*(*ndms)[i - 1] ) * ( *dmshifts )[nchans - 1] ) / ( tsamp ));
			*maxshift = (int) ceil((float) ( *maxshift + ( (float) ( SDIVINT*2*SNUMREG ) ) ) / (float) inBin[i - 1]) / (float) ( SDIVINT*2*SNUMREG );
			*maxshift = ( *maxshift ) * ( SDIVINT*2*SNUMREG ) * inBin[i - 1];
			if (( *maxshift ) > maxshift_high)
				maxshift_high = ( *maxshift );
		}
	}

	if (inBin[range - 1] > 1) {
		*maxshift = (int) ceil(( ( ( *dm_low )[range - 1] + ( *dm_step )[range - 1] * ( *ndms )[range - 1] ) * ( *dmshifts )[nchans - 1] ) / ( tsamp ));
		*maxshift = (int) ceil((float) ( *maxshift + ( (float) ( SDIVINT*2*SNUMREG ) ) ) / (float) inBin[range - 1]) / (float) ( SDIVINT*2*SNUMREG );
		*maxshift = *maxshift * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
		if (( *maxshift ) > maxshift_high)
			maxshift_high = ( *maxshift );
	}

	if (maxshift_high == 0)	{
		maxshift_high = (int) ceil(( ( ( *dm_low )[range - 1] + ( *dm_step )[range - 1] * ( ( *ndms )[range - 1] ) ) * ( *dmshifts )[nchans - 1] ) / tsamp);
	}
	*max_dm = ceil(( *dm_high )[range - 1]);

	*maxshift = ( maxshift_high + ( SNUMREG * 2 * SDIVINT ) );
	printf("\nRange:\t%d, MAXSHIFT:\t%d, Scrunch value:\t%d", range - 1, *maxshift, inBin[range - 1]);
	printf("\nMaximum dispersive delay:\t%.2f (s)", *maxshift * tsamp);

	if (*maxshift >= nsamp)	{
		printf("\n\nERROR!! Your maximum DM trial exceeds the number of samples you have.\nReduce your maximum DM trial\n\n");
		exit(1);
	}

	printf("\nDiagonal DM:\t%f", ( tsamp * nchans * 0.0001205 * powf(( fch1 + ( foff * ( nchans / 2 ) ) ), 3.0) ) / ( -foff * nchans ));
	if (*maxshift >= nsamp)	{
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
	( *t_processed ) = (int **) malloc(range * sizeof(int *));

	if (nchans < ( *max_ndms )) {
		// This means that we can cornerturn into the allocated output buffer 
		// without increasing the memory needed

		// Maximum number of samples we can fit in our GPU RAM is then given by:
		//max_tsamps = (unsigned int) ( (*gpu_memory) / ( sizeof(unsigned short) * ( (*max_ndms) + nchans ) ) ); // maximum number of timesamples we can fit into GPU memory
		max_tsamps = (unsigned int) ( (*gpu_memory) / ( sizeof(unsigned short)*nchans + sizeof(float)*(*max_ndms) + (size_t)(SPS_mem_requirement*MIN_DMS_PER_SPS_RUN ))); // maximum number of timesamples we can fit into GPU memory
		
		// Check that we dont have an out of range maxshift:
		if (( *maxshift ) > max_tsamps)	{
			printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
			printf("\nReduce your maximum dm or increase the size of your dm step");
			exit(0);
		}

		// Next check to see if nsamp fits in GPU RAM:
		if (nsamp < max_tsamps)	{
			// We have case 1)
			// Allocate memory to hold the values of nsamps to be processed
			unsigned long int local_t_processed = (unsigned long int) floor(( (double) ( nsamp - (*maxshift) ) / (double) inBin[range - 1] ) / (double) ( SDIVINT*2*SNUMREG )); //number of timesamples per block
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
			for (i = 0; i < range; i++)	{
				( *t_processed )[i] = (int *) malloc(sizeof(int)); // TODO: change to size_t
				( *t_processed )[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				( *t_processed )[i][0] = ( *t_processed )[i][0] * ( SDIVINT*2*SNUMREG );
			}
			( *num_tchunks ) = 1;
			printf("\nIn 1\n");
		}
		else {
			// We have case 3)
			// Work out how many time samples we can fit into ram 
			int samp_block_size = max_tsamps - ( *maxshift );

			// Work out how many blocks of time samples we need to complete the processing
			// upto nsamp-maxshift
			//int num_blocks = (int) floor(( (float) nsamp - ( *maxshift ) )) / ( (float) ( samp_block_size ) ) + 1;

			// Find the common integer amount of samples between all bins
			int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) inBin[range - 1] ) / (float) ( SDIVINT*2*SNUMREG ));
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
			
			int num_blocks = (int) floor(( (float) nsamp - (float)( *maxshift ) )) / ( (float) ( local_t_processed ) );

			// Work out the remaining fraction to be processed
			int remainder =  nsamp -  (num_blocks*local_t_processed ) - (*maxshift) ;
			remainder = (int) floor((float) remainder / (float) inBin[range - 1]) / (float) ( SDIVINT*2*SNUMREG );
			remainder = remainder * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];

			for (i = 0; i < range; i++)	{
				// Allocate memory to hold the values of nsamps to be processed
				( *t_processed )[i] = (int *) malloc((num_blocks + 1) * sizeof(int));
				// Remember the last block holds less!
				for (j = 0; j < num_blocks; j++) {
					( *t_processed )[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
					( *t_processed )[i][j] = ( *t_processed )[i][j] * ( SDIVINT*2*SNUMREG );
				}
				// fractional bit
				( *t_processed )[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				( *t_processed )[i][num_blocks] = ( *t_processed )[i][num_blocks] * ( SDIVINT*2*SNUMREG );
			}
			( *num_tchunks ) = num_blocks + 1;
			printf("\nIn 3\n");
			printf("\nnum_blocks:\t%d", num_blocks);
		}
	}
	else {
		// This means that we cannot cornerturn into the allocated output buffer 
		// without increasing the memory needed. Set the output buffer to be as large as the input buffer:

		// Maximum number of samples we can fit in our GPU RAM is then given by:
		//max_tsamps = (unsigned int) ( ( *gpu_memory ) / ( nchans * ( sizeof(float) + 2 * sizeof(unsigned short) ) ) );
		max_tsamps = (unsigned int) ( ( *gpu_memory ) / ( nchans * ( sizeof(float) + sizeof(unsigned short) )+ SPS_mem_requirement*MIN_DMS_PER_SPS_RUN ));

		// Check that we dont have an out of range maxshift:
		if (( *maxshift ) > max_tsamps) {
			printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
			printf("\nReduce your maximum dm or increase the size of your dm step");
			exit(0);
		}

		// Next check to see if nsamp fits in GPU RAM:
		if (nsamp < max_tsamps) {
			// We have case 2)
			// Allocate memory to hold the values of nsamps to be processed
			int local_t_processed = (int) floor(( (float) ( nsamp - ( *maxshift ) ) / (float) inBin[range - 1] ) / (float) ( SDIVINT*2*SNUMREG ));
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
			for (i = 0; i < range; i++) {
				( *t_processed )[i] = (int *) malloc(sizeof(int));
				( *t_processed )[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				( *t_processed )[i][0] = ( *t_processed )[i][0] * ( SDIVINT*2*SNUMREG );
			}
			( *num_tchunks ) = 1;
			printf("\nIn 2\n");
		}
		else {
			// We have case 4)
			// Work out how many time samples we can fit into ram 
			int samp_block_size = max_tsamps - ( *maxshift );

			// Work out how many blocks of time samples we need to complete the processing
			// upto nsamp-maxshift
			//int num_blocks = (int) floor(( (float) nsamp - (float) ( *maxshift ) ) / ( (float) samp_block_size ));

			// Find the common integer amount of samples between all bins
			int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) inBin[range - 1] ) / (float) ( SDIVINT*2*SNUMREG ));
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
			
			// samp_block_size was not used to calculate remainder instead there is local_t_processed which might be different
			int num_blocks = (int) floor(( (float) nsamp - (float) ( *maxshift ) ) / ( (float) local_t_processed ));

			// Work out the remaining fraction to be processed
			int remainder = nsamp - ( num_blocks * local_t_processed ) - ( *maxshift );
			remainder = (int) floor((float) remainder / (float) inBin[range - 1]) / (float) ( SDIVINT*2*SNUMREG );
			remainder = remainder * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];

			for (i = 0; i < range; i++)	{
				// Allocate memory to hold the values of nsamps to be processed
				( *t_processed )[i] = (int *) malloc(( num_blocks + 1 ) * sizeof(int));
				// Remember the last block holds less!
				for (j = 0; j < num_blocks; j++) {
					( *t_processed )[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
					( *t_processed )[i][j] = ( *t_processed )[i][j] * ( SDIVINT*2*SNUMREG );
				}
				// fractional bit
				( *t_processed )[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				( *t_processed )[i][num_blocks] = ( *t_processed )[i][num_blocks] * ( SDIVINT*2*SNUMREG );
			}
			( *num_tchunks ) = num_blocks + 1;
			printf("\nIn 4\n");
		}
	}
	printf("\nMaxshift memory needed:\t%lu MB", nchans * ( *maxshift ) * sizeof(unsigned short) / 1024 / 1024);
	if (nchans < ( *max_ndms ))	{
		printf("\nOutput memory needed:\t%lu MB", ( *max_ndms ) * ( *maxshift ) * sizeof(float) / 1024 / 1024);
	}
	else {
		printf("\nOutput memory needed:\t%lu MB", nchans * ( *maxshift ) * sizeof(float) / 1024 / 1024);
	}

}
