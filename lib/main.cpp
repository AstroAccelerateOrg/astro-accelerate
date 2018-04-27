#include "headers/headers_mains.h"
#include "headers/host_main_function.h" // Added by Nassim.O
#include "headers/host_debug.h"
#include "headers/host_get_user_input.h"
#include "headers/host_get_file_data.h"
#include "headers/host_allocate_memory.h"
#include "headers/host_get_recorded_data.h"
#include "headers/params.h"
#include "headers/device_init.h"


int main(int argc, char* argv[]) {
	// AstroAccelerate configuration. Which modules to run, debug options
	AA_Parameters AA_params;
	// DDTR data description
	DDTR_InputData DDTR_data;
	// DDTR plan
	DDTR_Plan DDTR_plan;
	//FDAS parameters
	FDAS_Parameters FDAS_params;
	//Single pulse search parameters
	SPS_Parameters SPS_params;
	//Periodicity search parameters
	PRS_Parameters PRS_params;
	
	// Internal code variables
	// File pointer
	FILE *fp = NULL;
	// Available memory for AA
	size_t gpu_memory = 0;
	// Input data for DDTR
	unsigned short *input_buffer = NULL;
	// Output data for DDTR
	float ***output_buffer = NULL;
	// Input data for DDTR on the device
	unsigned short *d_input = NULL;
	// Output data for DDTR on the device
	float *d_output = NULL;

	// Timing parameters
	clock_t start_time = clock();

	// Users desired de-dispersion strategy. Pick up user defined values from the CLI.
	// This reads user defined de-dispersion strategy from file and store it into DDTR plan, configures what modules AA should perform and configure these modules.
	get_user_input(&fp, argc, argv, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);
	// Input -> DDTR plan, All module parameters
				
	if (AA_params.enable_debug == 1)
		debug(1, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);
		
	// Reads telescope parameters from the header of the input file and then counts the number of samples in the input data file.
	DDTR_InputData DDTR_data;
	get_file_data(&fp, &DDTR_data);
	// This reads parameters of the data. Input -> time/frequency data type
		
	if (AA_params.enable_debug == 1)
		debug(3, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);

	// Allocate memory on host for input
	//DDTR_data.Allocate_memory_for_input();
	allocate_memory_cpu_input(&input_buffer, DDTR_data.nsamp, DDTR_data.nchans);
	  
	if (AA_params.enable_debug == 1)
		debug(5, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);

	// Store the recorded telescope data contained in the input filterbank file in the allocated memory.
	get_recorded_data(&fp, DDTR_data.nsamp, DDTR_data.nchans, DDTR_data.nbits, &input_buffer);
	// Reads raw data from disk. Input -> time/frequency data type
	
	if (AA_params.enable_debug == 1)
		debug(7, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);
	
	init_gpu(&gpu_memory, CARD);
	if(AA_params.enable_debug == 1) 
		debug(2, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);

	stratagy(&DDTR_plan, gpu_memory, DDTR_data, AA_params);
	
	if(AA_params.enable_debug == 1) 
		debug(4, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);
	
	
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
