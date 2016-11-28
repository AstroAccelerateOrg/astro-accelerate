#include "AstroAccelerate/DedispersionPlan.h"

namespace ska {
namespace astroaccelerate {
namespace sps {

	DedispersionPlan::DedispersionPlan()
	{
		_in_bin 			= NULL;
		_out_bin 			= NULL;
		_maxshift 		= 0;
		_dm_low 			= NULL;
		_dm_high 			= NULL;
		_dm_step 			= NULL;
		_dmshifts 		= NULL;
		_ndms 				= NULL;
		_max_ndms 		= 0;
		_range 				= 0;
		_t_processed 	= NULL;
		_nbits 				= 0;
		_nifs 				= 0;
		_tstart 			= 0.0f;
		_tsamp 				= 0.0f;
		_nsamp 				= 0;
		_nsamples 		= 0;
		_max_samps 		= 0;
		_nchans 			= 0;
		_fch1 				= 0.0f;
		_foff 				= 0.0f;
	}

	DedispersionPlan::~DedispersionPlan()
	{
	}
	void DedispersionPlan::get_file_data(FILE **fp)
		{
			/*
			fpos_t file_loc;

		char *string = (char *) malloc(80 * sizeof(char));

		int nchar;
		int nbytes = sizeof(int);

		unsigned long int total_data;

		double temp;

		while (1)
		{
			strcpy(string, "ERROR");
			//if (fread(&nchar, sizeof(int), 1, *fp) != 1)
			//	fprintf(stderr, "Error while reading file\n");
			fread(&nchar, sizeof(int), 1, *fp);
			if (feof(*fp))
				exit(0);

			if (nchar > 1 && nchar < 80)
			{
				//if (fread(string, nchar, 1, *fp) != 1)
				//	fprintf(stderr, "Error while reading file\n");
				fread(string, nchar, 1, *fp);
				string[nchar] = '\0';
				// For debugging only
				printf("\n%d\t%s", nchar, string), fflush(stdout);
				nbytes += nchar;

				if (strcmp(string, "HEADER_END") == 0)
					break;

				if (strcmp(string, "tsamp") == 0)
				{
					if (fread(&temp, sizeof(double), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
					*tsamp = (float) temp;
				}
				else if (strcmp(string, "tstart") == 0)
				{
					if (fread(&temp, sizeof(double), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
					*tstart = (float) temp;
				}
				else if (strcmp(string, "fch1") == 0)
				{
					if (fread(&temp, sizeof(double), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
					*fch1 = (float) temp;
				}
				else if (strcmp(string, "foff") == 0)
				{
					if (fread(&temp, sizeof(double), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
					*foff = (float) temp;
				}
				else if (strcmp(string, "nchans") == 0)
				{
					if (fread(nchans, sizeof(int), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
				}
				else if (strcmp(string, "nifs") == 0)
				{
					if (fread(nifs, sizeof(int), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
				}
				else if (strcmp(string, "nbits") == 0)
				{
					if (fread(nbits, sizeof(int), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
				}
				else if (strcmp(string, "nsamples") == 0)
				{
					if (fread(nsamples, sizeof(int), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
				}
			}
		}

		// Check that we are working with one IF channel
		if (*nifs != 1)
		{
			printf("\nERROR!! Can only work with one IF channel!\n");
			exit(1);
		}

		fgetpos(*fp, &file_loc);

		if (( *nbits ) == 32)
		{
			// Allocate a tempory buffer to store a line of frequency data
			float *temp_buffer = (float *) malloc(( *nchans ) * sizeof(float));

			// Count how many time samples we have
			total_data = 0;
			while (!feof(*fp))
			{
				//if (fread(temp_buffer, sizeof(float), ( *nchans ), *fp) != ( *nchans ))
				//	fprintf(stderr, "Error while reading file\n");
				fread(temp_buffer, sizeof(float), ( *nchans ), *fp);
				total_data++;
			}
			*nsamp = total_data - 1;

			free(temp_buffer);
		}
		else if (( *nbits ) == 8)
		{
			// Allocate a tempory buffer to store a line of frequency data
			unsigned char *temp_buffer = (unsigned char *) malloc(( *nchans ) * sizeof(unsigned char));

			total_data = 0;
			while (!feof(*fp))
			{
				//if (fread(temp_buffer, sizeof(unsigned char), ( *nchans ), *fp) != ( *nchans ))
				//	fprintf(stderr, "Error while reading file\n");
				fread(temp_buffer, sizeof(unsigned char), ( *nchans ), *fp);
				total_data++;
			}
			*nsamp = total_data - 1;

			free(temp_buffer);
		}
		else
		{
			printf("\n\n======================= ERROR =======================\n");
			printf(" Currently this code only runs with 1, 8 and 32 bit data\n");
			printf("\n=====================================================\n");
		}

		// Move the file pointer back to the end of the header
		fsetpos(*fp, &file_loc);

			 */
		}
	void DedispersionPlan::make_strategy(UserInput &user_input)
	{
/*
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

		if (power != 2.0)
		{
			// Calculate time independent dm shifts
			for (c = 0; c < nchans; c++)
			{
				( *dmshifts )[c] = 4148.741601f * ( ( 1.0 / pow(( fch1 + ( foff * c ) ), power) ) - ( 1.0 / pow(fch1, power) ) );
			}
		}
		else
		{
			// Calculate time independent dm shifts
			for (c = 0; c < nchans; c++)
			{
				( *dmshifts )[c] = (float) ( 4148.741601f * ( ( 1.0 / pow((double) ( fch1 + ( foff * c ) ), power) ) - ( 1.0 / pow((double) fch1, power) ) ) );
			}
		}

		for (i = 0; i < range; i++)
		{
			modff(( ( (int) ( ( user_dm_high[i] - user_dm_low[i] ) / user_dm_step[i] ) + SDIVINDM ) / SDIVINDM ), &n);
			( *ndms )[i] = (int) ( (int) n * SDIVINDM );
			if (*max_ndms < ( *ndms )[i])
				*max_ndms = ( *ndms )[i];
			*total_ndms = *total_ndms + ( *ndms )[i];
		}
		printf("\nMaximum number of dm trials in any of the range steps:\t%d", *max_ndms);

		( *dm_low )[0] = user_dm_low[0];
		( *dm_high )[0] = ( *ndms )[0] * ( user_dm_step[0] );
		( *dm_step )[0] = user_dm_step[0];
		for (i = 1; i < range; i++)
		{
			( *dm_low )[i] = ( *dm_high )[i - 1];
			( *dm_high )[i] = ( *dm_low )[i] + ( *ndms )[i] * user_dm_step[i];
			( *dm_step )[i] = user_dm_step[i];

			if (inBin[i - 1] > 1)
			{
				*maxshift = (int) ceil(( ( ( *dm_low )[i - 1] + ( *dm_step )[i - 1] * ( *ndms )[i - 1] ) * ( *dmshifts )[nchans - 1] ) / ( tsamp ));
				*maxshift = (int) ceil((float) ( *maxshift + ( (float) ( SDIVINT * ( SNUMREG ) ) ) ) / (float) inBin[i - 1]) / (float) ( SDIVINT * ( SNUMREG ) );
				*maxshift = ( *maxshift ) * ( SDIVINT * ( SNUMREG ) ) * inBin[i - 1];
				if (( *maxshift ) > maxshift_high)
					maxshift_high = ( *maxshift );
			}
		}

		if (inBin[range - 1] > 1)
		{
			*maxshift = (int) ceil(( ( ( *dm_low )[range - 1] + ( *dm_step )[range - 1] * ( *ndms )[range - 1] ) * ( *dmshifts )[nchans - 1] ) / ( tsamp ));
			*maxshift = (int) ceil((float) ( *maxshift + ( (float) ( SDIVINT * ( SNUMREG ) ) ) ) / (float) inBin[range - 1]) / (float) ( SDIVINT * ( SNUMREG ) );
			*maxshift = *maxshift * ( SDIVINT * ( SNUMREG ) ) * inBin[range - 1];
			if (( *maxshift ) > maxshift_high)
				maxshift_high = ( *maxshift );
		}

		if (maxshift_high == 0)
		{
			maxshift_high = (int) ceil(( ( ( *dm_low )[range - 1] + ( *dm_step )[range - 1] * ( ( *ndms )[range - 1] ) ) * ( *dmshifts )[nchans - 1] ) / tsamp);
		}
		*max_dm = ceil(( *dm_high )[range - 1]);

		*maxshift = ( maxshift_high + 2 * ( SNUMREG * SDIVINT ) );
		printf("\nRange:\t%d, MAXSHIFT:\t%d, Scrunch value:\t%d", range - 1, *maxshift, inBin[range - 1]);
		printf("\nMaximum dispersive delay:\t%.2f (s)", *maxshift * tsamp);

		if (*maxshift >= nsamp)
		{
			printf("\n\nERROR!! Your maximum DM trial exceeds the number of samples you have.\nReduce your maximum DM trial\n\n");
			exit(1);
		}

		printf("\nDiagonal DM:\t%f", ( tsamp * nchans * 0.0001205 * powf(( fch1 + ( foff * ( nchans / 2 ) ) ), 3.0) ) / ( -foff * nchans ));
		if (*maxshift >= nsamp)
		{
			printf("ERROR!! Your maximum DM trial exceeds the number of samples you have.\nReduce your maximum DM trial");
			exit(1);
		}

		/* Four cases:
		 * 1) nchans < max_ndms & nsamp fits in GPU RAM
		 * 2) nchans > max_ndms & nsamp fits in GPU RAM
		 * 3) nchans < max_ndms & nsamp does not fit in GPU RAM
		 * 4) nchans > max_ndms & nsamp does not fit in GPU RAM
		 */
/*
		int max_tsamps;

		// Allocate memory to store the t_processed ranges:
		( *t_processed ) = (int **) malloc(range * sizeof(int *));

		if (nchans < ( *max_ndms ))
		{
			// This means that we can cornerturn into the allocated output buffer
			// without increasing the memory needed

			// Maximum number of samples we can fit in our GPU RAM is then given by:
			max_tsamps = (int) ( ( *gpu_memory ) / ( sizeof(unsigned short) * ( ( *max_ndms ) + nchans ) ) );

			// Check that we dont have an out of range maxshift:
			if (( *maxshift ) > max_tsamps)
			{
				printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
				printf("\nReduce your maximum dm or increase the size of your dm step");
				exit(0);
			}

			// Next check to see if nsamp fits in GPU RAM:
			if (nsamp < max_tsamps)
			{
				// We have case 1)
				// Allocate memory to hold the values of nsamps to be processed
				int local_t_processed = (int) floor(( (float) ( nsamp - ( *maxshift ) ) / (float) inBin[range - 1] ) / (float) ( SDIVINT * ( SNUMREG ) ));
				local_t_processed = local_t_processed * ( SDIVINT * ( SNUMREG ) ) * inBin[range - 1];
				for (i = 0; i < range; i++)
				{
					( *t_processed )[i] = (int *) malloc(sizeof(int));
					( *t_processed )[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
					( *t_processed )[i][0] = ( *t_processed )[i][0] * ( SDIVINT * ( SNUMREG ) );
				}
				( *num_tchunks ) = 1;
				printf("\nIn 1\n");
			}
			else
			{
				// We have case 3)
				// Work out how many time samples we can fit into ram
				int samp_block_size = max_tsamps - ( *maxshift );
				//int samp_block_size = max_tsamps;

				// Work out how many blocks of time samples we need to complete the processing
				// upto nsamp-maxshift
				int num_blocks = (int) floor(( (float) nsamp - ( *maxshift ) )) / ( (float) ( samp_block_size ) ) + 1;

				// Find the common integer amount of samples between all bins
				int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) inBin[range - 1] ) / (float) ( SDIVINT * ( SNUMREG ) ));
				local_t_processed = local_t_processed * ( SDIVINT * ( SNUMREG ) ) * inBin[range - 1];

				// Work out the remaining fraction to be processed
				int remainder = ( nsamp - ( ( num_blocks - 1 ) * local_t_processed ) - ( *maxshift ) );
				remainder = (int) floor((float) remainder / (float) inBin[range - 1]) / (float) ( SDIVINT * ( SNUMREG ) );
				remainder = remainder * ( SDIVINT * ( SNUMREG ) ) * inBin[range - 1];

				for (i = 0; i < range; i++)
				{
					// Allocate memory to hold the values of nsamps to be processed
					( *t_processed )[i] = (int *) malloc(num_blocks * sizeof(int));
					// Remember the last block holds less!
					for (j = 0; j < num_blocks - 1; j++)
					{
						( *t_processed )[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
						( *t_processed )[i][j] = ( *t_processed )[i][j] * ( SDIVINT * ( SNUMREG ) );
					}
					// fractional bit
					( *t_processed )[i][num_blocks - 1] = (int) floor(( (float) ( remainder ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
					( *t_processed )[i][num_blocks - 1] = ( *t_processed )[i][num_blocks - 1] * ( SDIVINT * ( SNUMREG ) );
				}
				( *num_tchunks ) = num_blocks;
				printf("\nIn 3\n");
				printf("\nnum_blocks:\t%d", num_blocks);
			}
		}
		else
		{
			// This means that we cannot cornerturn into the allocated output buffer
			// without increasing the memory needed. Set the output buffer to be as large as the input buffer:

			// Maximum number of samples we can fit in our GPU RAM is then given by:
			max_tsamps = (int) ( ( *gpu_memory ) / ( nchans * ( sizeof(float) + 2 * sizeof(unsigned short) ) ) );

			// Check that we dont have an out of range maxshift:
			if (( *maxshift ) > max_tsamps)
			{
				printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
				printf("\nReduce your maximum dm or increase the size of your dm step");
				exit(0);
			}

			// Next check to see if nsamp fits in GPU RAM:
			if (nsamp < max_tsamps)
			{
				// We have case 2)
				// Allocate memory to hold the values of nsamps to be processed
				int local_t_processed = (int) floor(( (float) ( nsamp - ( *maxshift ) ) / (float) inBin[range - 1] ) / (float) ( SDIVINT * ( SNUMREG ) ));
				local_t_processed = local_t_processed * ( SDIVINT * ( SNUMREG ) ) * inBin[range - 1];
				for (i = 0; i < range; i++)
				{
					( *t_processed )[i] = (int *) malloc(sizeof(int));
					( *t_processed )[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
					( *t_processed )[i][0] = ( *t_processed )[i][0] * ( SDIVINT * ( SNUMREG ) );
				}
				( *num_tchunks ) = 1;
				printf("\nIn 2\n");
			}
			else
			{
				// We have case 4)
				// Work out how many time samples we can fit into ram
				int samp_block_size = max_tsamps - ( *maxshift );

				// Work out how many blocks of time samples we need to complete the processing
				// upto nsamp-maxshift
				int num_blocks = (int) floor(( (float) nsamp - (float) ( *maxshift ) ) / ( (float) samp_block_size ));

				// Find the common integer amount of samples between all bins
				int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) inBin[range - 1] ) / (float) ( SDIVINT * ( SNUMREG ) ));
				local_t_processed = local_t_processed * ( SDIVINT * ( SNUMREG ) ) * inBin[range - 1];

				// Work out the remaining fraction to be processed
				int remainder = nsamp - ( num_blocks * local_t_processed ) - ( *maxshift );

				for (i = 0; i < range; i++)
				{
					// Allocate memory to hold the values of nsamps to be processed
					( *t_processed )[i] = (int *) malloc(( num_blocks + 1 ) * sizeof(int));
					// Remember the last block holds less!
					for (j = 0; j < num_blocks; j++)
					{
						( *t_processed )[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
						( *t_processed )[i][j] = ( *t_processed )[i][j] * ( SDIVINT * ( SNUMREG ) );
					}
					// fractional bit
					( *t_processed )[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
					( *t_processed )[i][num_blocks] = ( *t_processed )[i][num_blocks] * ( SDIVINT * ( SNUMREG ) );
				}
				( *num_tchunks ) = num_blocks + 1;
				printf("\nIn 4\n");
			}
		}

		printf("\nMaxshift memory needed:\t%lu MB", nchans * ( *maxshift ) * sizeof(unsigned short) / 1024 / 1024);
		if (nchans < ( *max_ndms ))
		{
			printf("\nOutput memory needed:\t%lu MB", ( *max_ndms ) * ( *maxshift ) * sizeof(float) / 1024 / 1024);
		}
		else
		{
			printf("\nOutput memory needed:\t%lu MB", nchans * ( *maxshift ) * sizeof(float) / 1024 / 1024);
		}
*/
	}


} // namespace sps
} // namespace astroaccelerate
} // namespace ska
