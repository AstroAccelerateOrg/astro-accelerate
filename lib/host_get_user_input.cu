/* This module recieves user input from the CLI and then opens the input data file.
 * The input data file should always be the last argument on the CLI.
 */

#include <stdio.h>
#include <wordexp.h>
#include "headers/params.h"
#include "headers/host_help.h"
#include "headers/headers_mains.h"

void get_user_input(FILE **fp, int argc, char *argv[], DDTR_Plan *DDTR_plan, AA_Parameters *AA_params, MSD_Parameters *MSD_params, SPS_Parameters *SPS_params, PRS_Parameters *PRS_params, FDAS_Parameters *FDAS_params) {

	FILE *fp_in = NULL;
	int nRanges, error;
	//int nb_selected_dm;

	char string[100];
	int i;

	error = 0;
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
		nRanges = 0;
		//nb_selected_dm = 0;
		while (!feof(fp_in))
		{
			if ( fscanf(fp_in, "%s", string) == 0 )
			{
				fprintf(stderr, "failed to read input\n");
				exit(0);
			}
			if (strcmp(string, "range") == 0)
				nRanges++;
			// Ignored as it is not used
			//if (strcmp(string, "selected_dm") == 0)
			//	( *nb_selected_dm )++;
		}
		rewind(fp_in);
		
		DDTR_plan->nRanges = nRanges;
		error = DDTR_plan->Allocate_user_ranges();
		if(error>1) exit(98);
		
		//*user_dm_low = (float *) malloc(( *range ) * sizeof(float));
		//*user_dm_high = (float *) malloc(( *range ) * sizeof(float));
		//*user_dm_step = (float *) malloc(( *range ) * sizeof(float));
		//*inBin = (int *) malloc(( *range ) * sizeof(int));

		// temporary variables to read dm range
		float temp_low  = 0;
		float temp_high = 0;
		float temp_step = 0;
		int temp_in_bin = 0;
		int temp_out_bin= 0;

		// read dm range if enabled
		i=0;
		while (!feof(fp_in))
		{
			if ( fscanf(fp_in, "%s %f %f %f %d %d\n", string, &temp_low, &temp_high, &temp_step, &temp_in_bin, &temp_out_bin) == 0 )
			{
				fprintf(stderr, "failed to read input\n");
				exit(0);
			}
			if (strcmp(string, "range") == 0) {
				DDTR_plan->user_dm_low[i]  = temp_low;
				DDTR_plan->user_dm_high[i] = temp_high;
				DDTR_plan->user_dm_step[i] = temp_step;
				DDTR_plan->inBin[i] = temp_in_bin;
				i++;
			}
		}
		rewind(fp_in);

		//
		// Ignored
		//*selected_dm_low  = (float *) malloc(( *nb_selected_dm ) * sizeof(float));
		//*selected_dm_high = (float *) malloc(( *nb_selected_dm ) * sizeof(float));
		//// temporary variables to read selected dm range
		//float temp_selected_low  = 0;
		//float temp_selected_high = 0;
		//// read selected dm range if enabled
		//i=0;
		//while (!feof(fp_in))
		//{
		//	if ( fscanf(fp_in, "%s %f %f\n", string, &temp_selected_low, &temp_selected_high) == 0 )
		//	{
		//		fprintf(stderr, "failed to read input\n");
		//		exit(0);
		//	}
		//	if (strcmp(string, "selected_dm") == 0)
		//	{
		//		(*selected_dm_low)[i]  = temp_selected_low;
		//		(*selected_dm_high)[i] = temp_selected_high;
		//		i++;
		//	}
		//}
		//rewind(fp_in);

		while (!feof(fp_in))
		{
			if ( fscanf(fp_in, "%s", string) == 0 )
			{
				fprintf(stderr, "failed to read input\n");
				exit(0);
			}
			if (strcmp(string, "debug") == 0)
				AA_params->enable_debug = 1;
			if (strcmp(string, "analysis") == 0)
				AA_params->enable_analysis = 1;
			if (strcmp(string, "periodicity") == 0)
				AA_params->enable_periodicity = 1;
			if (strcmp(string, "acceleration") == 0)
				AA_params->enable_acceleration = 1;
			if (strcmp(string, "zero_dm") == 0)
				AA_params->enable_zero_dm = 1;
			if (strcmp(string, "zero_dm_with_outliers") == 0)
				AA_params->enable_zero_dm_with_outliers = 1;
			if (strcmp(string, "rfi") == 0)
				AA_params->enable_rfi = 1;
			if (strcmp(string, "oldrfi") == 0)
				AA_params->enable_old_rfi = 1;
			if (strcmp(string, "analysis_debug") == 0)
				AA_params->analysis_debug = 1;
			if (strcmp(string, "failsafe") == 0)
				AA_params->failsafe = 1;
			
			if (strcmp(string, "threshold") == 0) {
				FDAS_params->candidate_algorithm = 1;
				SPS_params->candidate_algorithm = 1;
				PRS_params->candidate_algorithm = 1;
			}
			
			
			if (strcmp(string, "output_ffdot_plan") == 0)
				FDAS_params->enable_output_ffdot_plan = 1;
			if (strcmp(string, "output_fdas_list") == 0)
				FDAS_params->enable_output_fdas_list = 1;
			if (strcmp(string, "fdas_custom_fft") == 0)
				FDAS_params->enable_fdas_custom_fft = 1;
			if (strcmp(string, "fdas_inbin") == 0)
				FDAS_params->enable_fdas_inbin = 1;
			if (strcmp(string, "fdas_norm") == 0)
				FDAS_params->enable_fdas_norm = 1;

			
			if (strcmp(string, "baselinenoise") == 0)
				MSD_params->enable_outlier_rejection = 1;
			if (strcmp(string, "sigma_constant") == 0) {
				if ( fscanf(fp_in, "%f", &(MSD_params->OR_sigma_multiplier)) == 0 ) {
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			

			//if (strcmp(string, "output_dmt") == 0)
			//	*output_dmt = 1;			
			//if (strcmp(string, "multi_file") == 0)
			//	*multi_file = 1;
			if (strcmp(string, "max_boxcar_width_in_sec") == 0) {
				if ( fscanf(fp_in, "%f", &(SPS_params->max_boxcar_width_in_sec)) == 0 ) {
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			if (strcmp(string, "sigma_cutoff") == 0) {
				float ftemp = 0;
				if ( fscanf(fp_in, "%f", &ftemp) == 0) {
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
				SPS_params->sigma_cutoff = ftemp;
				FDAS_params->sigma_cutoff = ftemp;
			}
			
			//-------------------> Periodicity search parameters
			if (strcmp(string, "periodicity_sigma_cutoff") == 0) {
				if ( fscanf(fp_in, "%f", &(PRS_params->sigma_cutoff)) == 0 ) {
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			if (strcmp(string, "periodicity_harmonics") == 0) {
				if ( fscanf(fp_in, "%d", &(PRS_params->nHarmonics)) == 0 ) {
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			//-------------------<
			
			if (strcmp(string, "power") == 0) {
				if ( fscanf(fp_in, "%f", &(DDTR_plan->power)) == 0 ) {
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			
			/*
			if (strcmp(string, "narrow") == 0)
			{
				if ( fscanf(fp_in, "%f", narrow) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			if (strcmp(string, "wide") == 0)
			{
				if ( fscanf(fp_in, "%f", wide) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			if (strcmp(string, "nboots") == 0)
			{
				if ( fscanf(fp_in, "%d", nboots) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			if (strcmp(string, "navdms") == 0)
			{
				if ( fscanf(fp_in, "%d", navdms) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			if (strcmp(string, "nwindows") == 0)
			{
				if ( fscanf(fp_in, "%d", ntrial_bins) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			if (strcmp(string, "nsearch") == 0)
			{
				if ( fscanf(fp_in, "%d", nsearch) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			if (strcmp(string, "aggression") == 0)
			{
				if ( fscanf(fp_in, "%f", aggression) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}
			*/
			
			if (strcmp(string, "file") == 0) {
				// this command expand "~" to "home/username/"
			    wordexp_t expanded_string;

				if ( fscanf(fp_in, "%s", string) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
				wordexp(string, &expanded_string, 0);
			    if (( *fp = fopen(expanded_string.we_wordv[0], "rb") ) == NULL)
				{
					fprintf(stderr, "Invalid data file!\n");
					exit(0);
				}
				wordfree(&expanded_string);
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
