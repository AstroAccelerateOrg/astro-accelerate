/* This module recieves user input from the CLI and then opens the input data file.
 * The input data file should always be the last argument on the CLI.
 */

#include <stdio.h>
#include "AstroAccelerate/params.h"
#include "AstroAccelerate/host_help.h"

void get_user_input(FILE **fp, int argc, char *argv[], int *multi_file, int *enable_debug, int *enable_analysis, int *enable_periodicity, int *enable_acceleration, int *output_dmt, int *enable_zero_dm, int *nboots, int *ntrial_bins, int *navdms, float *narrow, float *wide, float *aggression, int *nsearch, int **inBin, int **outBin, float *power, float *sigma_cutoff, int *range, float **user_dm_low, float **user_dm_high, float **user_dm_step)
{

	FILE *fp_in = NULL;

	char string[100];
	int i;

	//{{{ Read in the command line parameters and open the input file

	if (argc < 2)
	{
		fprintf(stderr, "Need input file.\n");
		exit(0);
	}
	else if (argc == 2 && strcmp(argv[1], "-help") != 0)
	{
		if (( fp_in = fopen(argv[1], "r") ) == NULL)
		{
			fprintf(stderr, "Invalid input file!\n");
			exit(0);
		}
		( *range ) = 0;
		while (!feof(fp_in))
		{
			fscanf(fp_in, "%s", string);
			if (strcmp(string, "range") == 0)
				( *range )++;
		}
		rewind(fp_in);

		*user_dm_low = (float *) malloc(( *range ) * sizeof(float));
		*user_dm_high = (float *) malloc(( *range ) * sizeof(float));
		*user_dm_step = (float *) malloc(( *range ) * sizeof(float));
		*outBin = (int *) malloc(( *range ) * sizeof(int));
		*inBin = (int *) malloc(( *range ) * sizeof(int));

		for (i = 0; i < *range; i++)
		{
			fscanf(fp_in, "%s %f %f %f %d %d\n", string, &( *user_dm_low )[i], &( *user_dm_high )[i], &( *user_dm_step )[i], &( *inBin )[i], &( *outBin )[i]);
		}

		rewind(fp_in);
		while (!feof(fp_in))
		{
			fscanf(fp_in, "%s", string);
			if (strcmp(string, "debug") == 0)
				*enable_debug = 1;
			if (strcmp(string, "analysis") == 0)
				*enable_analysis = 1;
			if (strcmp(string, "periodicity") == 0)
				*enable_periodicity = 1;
			if (strcmp(string, "acceleration") == 0)
				*enable_acceleration = 1;
			if (strcmp(string, "output_dmt") == 0)
				*output_dmt = 1;
			if (strcmp(string, "zero_dm") == 0)
				*enable_zero_dm = 1;
			if (strcmp(string, "multi_file") == 0)
				*multi_file = 1;
			if (strcmp(string, "sigma_cutoff") == 0)
				fscanf(fp_in, "%f", sigma_cutoff);
			if (strcmp(string, "narrow") == 0)
				fscanf(fp_in, "%f", narrow);
			if (strcmp(string, "wide") == 0)
				fscanf(fp_in, "%f", wide);
			if (strcmp(string, "nboots") == 0)
				fscanf(fp_in, "%d", nboots);
			if (strcmp(string, "navdms") == 0)
				fscanf(fp_in, "%d", navdms);
			if (strcmp(string, "nwindows") == 0)
				fscanf(fp_in, "%d", ntrial_bins);
			if (strcmp(string, "nsearch") == 0)
				fscanf(fp_in, "%d", nsearch);
			if (strcmp(string, "aggression") == 0)
				fscanf(fp_in, "%f", aggression);
			if (strcmp(string, "power") == 0)
				fscanf(fp_in, "%f", power);
			if (strcmp(string, "file") == 0)
			{
				fscanf(fp_in, "%s", string);
				if (( *fp = fopen(string, "rb") ) == NULL)
				{
					fprintf(stderr, "Invalid data file!\n");
					exit(0);
				}
			}
		}

	}
	else if (argc == 2 && strcmp(argv[1], "-help") == 0)
	{
		help();
	}
	else
	{
		fprintf(stderr, "Cannot recognise input, try \"./astro-accelerate -help.\"\n");
		exit(0);
	}
	//}}}
}
