#ifndef __GETUSERINPUT__
#define __GETUSERINPUT__

void get_user_input(FILE **fp, int argc, char *argv[], int *multi_file, int *enable_debug, int *enable_analysis, int *enable_periodicity, int *enable_acceleration, int *output_dmt, int *enable_zero_dm, int *nboots, int *ntrial_bins, int *navdms, float *narrow, float *wide, float *aggression, int *nsearch, int **inBin, int **outBin, float *power, float *sigma_cutoff, int *range, float **user_dm_low, float **user_dm_high, float **user_dm_step);
#endif
