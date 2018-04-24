#include "headers/headers_mains.h"
#include "headers/host_main_function.h" // Added by Nassim.O
#include "headers/host_debug.h"
#include "headers/host_get_user_input.h"
#include "headers/host_get_file_data.h"
#include "headers/host_allocate_memory.h"
#include "headers/host_get_recorded_data.h"
#include "headers/params.h"

int main(int argc, char* argv[])
{
	// Internal code variables
	// File pointers
	FILE *fp = NULL;
	// Counters and flags
	int i, t, dm_range;
	int range = 0;
	int nb_selected_dm = 0;
	int enable_debug = 0;
	int enable_analysis = 0;
	int enable_acceleration = 0;
	int enable_periodicity = 0;
	int output_dmt = 0;
	int enable_zero_dm = 0;
	int enable_zero_dm_with_outliers = 0;
	int enable_rfi = 0;
	int enable_old_rfi = 0;
	int enable_sps_baselinenoise=0;
	int enable_fdas_custom_fft = 0;
	int enable_fdas_inbin = 0;
	int enable_fdas_norm = 0;
	int enable_output_ffdot_plan = 0;
	int enable_output_fdas_list = 0;
	int analysis_debug = 0;
    int *inBin = NULL;
	int *outBin = NULL;
	int *ndms = NULL;
	int maxshift = 0;
	int max_ndms = 0;
	int max_samps = 0;
	int num_tchunks = 0;
	int total_ndms = 0;
	int multi_file = 1;
	float max_dm = 0.0f;
	int candidate_algorithm=0;
	int failsafe = 0;
	// Memory sizes and pointers
	size_t inputsize = 0;
	size_t outputsize = 0;
	size_t gpu_inputsize = 0;
	size_t gpu_outputsize = 0;
	size_t gpu_memory = 0;
	unsigned short *input_buffer = NULL;
	float ***output_buffer = NULL;
	unsigned short *d_input = NULL;
	float *d_output = NULL;
	float *dmshifts = NULL;
	float *user_dm_low = NULL;
	float *user_dm_high = NULL;
	float *user_dm_step = NULL;
	float *dm_low = NULL;
	float *dm_high = NULL;
	float *dm_step = NULL;
	float *selected_dm_low = NULL;
	float *selected_dm_high = NULL;
	// Telescope parameters
	int nchans = 0;
	int nsamp = 0;
	int nbits = 0;
	int nsamples = 0;
	int nifs = 0;
	int **t_processed;
	int nboots = -1;
	int ntrial_bins;
	int navdms = 1;
	int nsearch = 3;
	float aggression = 2.5;
	float narrow = 0.001f;
	float wide = 0.1f;
	int maxshift_original;
	double tsamp_original;
	long int inc = 0;
	float tstart = 0.0f;
	float tstart_local = 0.0f;
	float tsamp = 0.0f;
	float fch1 = 0.0f;
	float foff = 0.0f;
	// Analysis variables
	float power = 2.0f;
	float sigma_cutoff = 6.0f;
	float sigma_constant = 4.0f;
	float max_boxcar_width_in_sec = 0.5f;
	// Periodicity search
	float periodicity_sigma_cutoff = 6;
	int periodicity_nHarmonics = 32;

	// Timing parameters
	clock_t start_time = clock();

	// Users desired de-dispersion strategy. Pick up user defined values from the CLI.
	get_user_input(&fp, argc, argv, &multi_file, &enable_debug, &enable_analysis, &enable_periodicity, &enable_acceleration, &enable_output_ffdot_plan, &enable_output_fdas_list, &output_dmt, &enable_zero_dm, &enable_zero_dm_with_outliers, &enable_rfi, &enable_old_rfi, &enable_fdas_custom_fft, &enable_fdas_inbin, &enable_fdas_norm, &nboots, &ntrial_bins, &navdms, &narrow, &wide, &aggression, &nsearch, &inBin, &outBin, &power, &sigma_cutoff, &sigma_constant, &max_boxcar_width_in_sec, &range, &user_dm_low, &user_dm_high, &user_dm_step, &candidate_algorithm, &enable_sps_baselinenoise, &selected_dm_low, &selected_dm_high, &nb_selected_dm, &analysis_debug, &failsafe, &periodicity_sigma_cutoff, &periodicity_nHarmonics);
	// This reads DDTR plan, configures what modules AA should perform and configure these modules.
	// Input -> DDTR plan, All module parameters
				
	if (AA_params.enable_debug == 1)
		debug(1, DDTR_data, DDTR_plan, AA_params, MSD_params, SPS_params, PRS_params, FDAS_params);
		
	// Reads telescope parameters from the header of the input file and then counts the number of samples in the input data file.
	DDTR_Data DDTR_data;
	get_file_data(&fp, &DDTR_data);
	// This reads parameters of the data. Input -> time/frequency data type
		
	if (AA_params.enable_debug == 1)
		debug(3, DDTR_data, DDTR_plan, AA_params, MSD_params, SPS_params, PRS_params, FDAS_params);

	// Allocate memory on host.
	allocate_memory_cpu_input(&fp, gpu_memory, maxshift, num_tchunks, max_ndms,
	  total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer,
	  &output_buffer, &d_input, &d_output, &gpu_inputsize, &gpu_outputsize,
	  &inputsize, &outputsize);
	  
	if (AA_params.enable_debug == 1)
		debug(5, DDTR_data, DDTR_plan, AA_params, MSD_params, SPS_params, PRS_params, FDAS_params);

	// Store the recorded telescope data contained in the input filterbank file in the allocated memory.
	get_recorded_data(&fp, nsamp, nchans, nbits, &input_buffer, &inputsize);
	// Reads raw data from disk. Input -> time/frequency data type
	
	if (AA_params.enable_debug == 1)
		debug(7, DDTR_data, DDTR_plan, AA_params, MSD_params, SPS_params, PRS_params, FDAS_params);
	
	// Check GPU available
	// Run stratagy
	// Allocate host output
		

	main_function (
	  argc, argv,
	  // Internal code variables
	  // File pointers
	  fp,
	  // Counters and flags
	  i, t, dm_range, range, enable_debug, enable_analysis, enable_acceleration, enable_output_ffdot_plan,
	  enable_output_fdas_list, enable_periodicity, output_dmt, enable_zero_dm, enable_zero_dm_with_outliers,
	  enable_rfi, enable_old_rfi, enable_sps_baselinenoise, enable_fdas_custom_fft, enable_fdas_inbin, enable_fdas_norm, inBin,
	  outBin, ndms, maxshift, max_ndms, max_samps, num_tchunks, total_ndms, multi_file, max_dm,
	  // Memory sizes and pointers
	  inputsize, outputsize, gpu_inputsize, gpu_outputsize, gpu_memory,
	  input_buffer, output_buffer, d_input, d_output, dmshifts, user_dm_low,
	  user_dm_high, user_dm_step, dm_low, dm_high, dm_step,
	  // Telescope parameters
	  nchans, nsamp, nbits, nsamples, nifs, t_processed, nboots, ntrial_bins,
	  navdms, nsearch, aggression, narrow, wide, maxshift_original,
	  tsamp_original, inc, tstart, tstart_local, tsamp, fch1, foff,
	  // Analysis variables
	  power, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, start_time, candidate_algorithm,
	  nb_selected_dm, selected_dm_low, selected_dm_high, analysis_debug, failsafe,
	  // Periodicity search
	  periodicity_sigma_cutoff, periodicity_nHarmonics
	);

	// write output here, not in the library

	fclose(fp);
	


	/*	
	template<typename ValueType>
	DmTime<ValueType>::~DmTime() {
		for(int i = 0; i < _range; ++i)
		{
			for(int j = 0; j < _ndms[i]; ++j)
			{
				free(_data[i][j]);
			}
			free(_data[i]);
		}
		free(_data);
	}
	*/

	free(output_buffer);
	free(t_processed);
	free(dm_low);
	free(dm_high);
	free(dm_step);
	free(dmshifts);
	free(user_dm_low);
	free(user_dm_high);
	free(user_dm_step);
	if (nb_selected_dm > 0)
	{
		free(selected_dm_low);
		free(selected_dm_high);
	}

	return 0;

}
