#include "headers/headers_mains.h"
#include "headers/host_main_function.h" // Added by Nassim.O
#include "headers/host_debug.h"
#include "headers/host_get_user_input.h"
#include "headers/host_get_file_data.h"
#include "headers/host_allocate_memory.h"
#include "headers/host_get_recorded_data.h"
#include "headers/host_stratagy.h"
#include "headers/params.h"
#include "headers/device_init.h"


int main(int argc, char* argv[]) {
	// Configurations
	// AstroAccelerate Which modules to run, debug options
	AA_Parameters AA_params;
	// DDTR plan
	DDTR_Plan DDTR_plan;
	//FDAS parameters
	FDAS_Parameters FDAS_params;
	//Single pulse search parameters
	SPS_Parameters SPS_params;
	//Periodicity search parameters
	PRS_Parameters PRS_params;
	//MSD parameters
	MSD_Parameters MSD_params;
	
	// DDTR input data description
	DDTR_InputData DDTR_data;
	// Input data for DDTR
	unsigned short *input_buffer = NULL;
	
	// DDTR output data. Description so far through DDTR_Plan.
	float ***output_buffer = NULL;
	// SPS candidate list
	float *h_SPS_candidatelist;
	// SPS number of candidates
	size_t nSPScandidates;
	
	// Internal code variables
	// File pointer
	FILE *fp = NULL;
	// Available memory for AA
	size_t gpu_memory = 0;
	
	// Timing parameters
	clock_t start_time = clock();

	// Users desired de-dispersion strategy. Pick up user defined values from the CLI.
	// This reads user defined de-dispersion strategy from file and store it into DDTR plan, configures what modules AA should perform and configure these modules.
	get_user_input(&fp, argc, argv, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);
	// Input -> DDTR plan, All module parameters
				
	if (AA_params.enable_debug == 1)
		debug(1, start_time, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);
		
	// Reads telescope parameters from the header of the input file and then counts the number of samples in the input data file.
	get_file_data(&fp, &DDTR_data);
	// This reads parameters of the data. Input -> time/frequency data type
		
	if (AA_params.enable_debug == 1)
		debug(3, start_time, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);

	// Allocate memory on host for input
	//DDTR_data.Allocate_memory_for_input();
	allocate_memory_cpu_input(&input_buffer, DDTR_data.nsamp, DDTR_data.nchans);
	  
	//if (AA_params.enable_debug == 1)
	//	debug(5, start_time, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);

	// Store the recorded telescope data contained in the input filterbank file in the allocated memory.
	get_recorded_data(&fp, DDTR_data.nsamp, DDTR_data.nchans, DDTR_data.nbits, &input_buffer);
	fclose(fp);
	
	if (AA_params.enable_debug == 1)
		debug(7, start_time, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);
	
	init_gpu(&gpu_memory, CARD);
	if(AA_params.enable_debug == 1) 
		debug(2, start_time, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);

	stratagy(&DDTR_plan, gpu_memory, &DDTR_data, &AA_params);
	
	if(AA_params.enable_debug == 1) 
		debug(4, start_time, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);
	
	// Allocate memory on host and device.
	allocate_memory_cpu_output(&output_buffer, &DDTR_plan);
	
	if(AA_params.enable_debug == 1) 
		debug(5, start_time, &DDTR_data, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params);
		
	
	main_function(&h_SPS_candidatelist, &nSPScandidates, output_buffer, input_buffer, &DDTR_plan, &AA_params, &MSD_params, &SPS_params, &PRS_params, &FDAS_params, start_time);
	
	
	// deallocate host output
	for(int i = 0; i < DDTR_plan.nRanges; i++) {
		for(int j = 0; j < DDTR_plan.ndms[i]; j++) {
			free(output_buffer[i][j]);
		}
		free(output_buffer[i]);
	}
	free(output_buffer);
	free(input_buffer);

	return 0;

}
