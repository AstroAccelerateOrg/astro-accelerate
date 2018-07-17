#include "headers/device_bin.h"
#include "headers/device_init.h"
#include "headers/device_dedisperse.h"
#include "headers/device_dedispersion_kernel.h"
#include "headers/device_zero_dm.h"
#include "headers/device_zero_dm_outliers.h"
#include "headers/device_rfi.h"
#include "headers/device_load_data.h"
#include "headers/device_corner_turn.h"

#include "headers/headers_mains.h"
#include "headers/device_DDTR_DedispersionTransform.h"

#include <helper_cuda.h>

#include "headers/device_AA_Parameters.h"


// MSD
#include "headers/device_MSD_Parameters.h"
#include "headers/device_MSD_Configuration.h"
#include "headers/device_MSD.h"
#include "headers/device_MSD_plane_profile.h"

// SPS
#include "headers/device_SPS_Search.h"
#include "headers/device_SPS_Search_AA.h"
#include "headers/device_SPS_DataDescription.h"
#include "headers/device_SPS_Parameters.h"
#include "headers/device_SPS_inplace_kernel.h" //Added by KA
#include "headers/device_SPS_inplace.h" //Added by KA
#include "headers/device_SNR_limited.h" //Added by KA
#include "headers/device_SPS_long.h" //Added by KA
#include "headers/device_threshold.h" //Added by KA
#include "headers/device_single_FIR.h" //Added by KA
#include "headers/device_analysis.h" //Added by KA
#include "headers/device_periods.h" //Added by KA
#include "headers/device_peak_find.h" //Added by KA

// PRS
#include "headers/device_power.h"
#include "headers/device_harmonic_summing.h"

#include "headers/device_save_data.h"
#include "headers/host_acceleration.h"
#include "headers/host_allocate_memory.h"
#include "headers/host_analysis.h"
#include "headers/host_export.h"
#include "headers/host_periods.h"
#include "headers/host_debug.h"
#include "headers/host_get_file_data.h"
#include "headers/host_get_recorded_data.h"
#include "headers/host_get_user_input.h"
#include "headers/host_help.h"
#include "headers/host_rfi.h"
#include "headers/host_stratagy.h"
#include "headers/host_write_file.h"

// fdas
#include "headers/device_acceleration_fdas.h"

#include "headers/host_main_function.h"

#include "headers/params.h"

#include "timer.h"


//#define EXPORT_DD_DATA

void main_function (
		SPS_CandidateList *candidatelist,
		float ***h_output_buffer, //---> AA Output data on the HOST
		unsigned short  *h_input_buffer, //---> AA Input data on the HOST
		DDTR_Plan *DDTR_plan,
		AA_Parameters *AA_params,
		MSD_Parameters *MSD_params,
		SPS_Parameters *SPS_params,
		PRS_Parameters *PRS_params,
		FDAS_Parameters *FDAS_params,
		clock_t start_time // Time measurements
		) {
	
	DDTR_DedispersionTransform DDTR;
	
	unsigned short *d_DDTR_input;
	float *d_DDTR_output;
	
	// Allocate memory on host and device. This should be inside some DDTR class not exposed to user.
	allocate_memory_gpu(&d_DDTR_input, &d_DDTR_output, DDTR_plan);
	
	if(AA_params->enable_debug == 1) {
		DDTR_InputData DDTR_data;
		debug(5, start_time, &DDTR_data, DDTR_plan, AA_params, MSD_params, SPS_params, PRS_params, FDAS_params);
	}

	checkCudaErrors(cudaGetLastError());
	
	// Clip RFI
	// This is done before transfering data to device (loading to DDTR) so no need to include it into DDTR class
	if (AA_params->enable_rfi) {
		printf("\nPerforming new CPU rfi...");
		rfi(DDTR_plan->nsamp, DDTR_plan->nchans, &h_input_buffer);
	}
	
	printf("\nDe-dispersing...");
	GpuTimer timer;
	timer.Start();


	//float tsamp_original = DDTR_plan->tsamp;
	//float local_tsamp = DDTR_plan->tsamp;
	//int maxshift_original = DDTR_plan->maxshift;
	//size_t local_maxshift = DDTR_plan->maxshift;
	
	// Need to keep these two
	float tstart_local = 0;
	size_t processed_samples = 0; //inc = number of processed time samples so far.
	

	while( !DDTR.endOfDDTRPlan(DDTR_plan) ){
		int DDTR_status = 0;
		
		// De-dispersion transform
		DDTR_status = DDTR.dedisperseNextChunk(d_DDTR_output, d_DDTR_input, h_input_buffer, DDTR_plan, AA_params);
		processed_samples = DDTR.getProcessedSamples();
		tstart_local = ( DDTR_plan->tsamp*processed_samples );
		
		if(DDTR_status<0) printf("Damn! This should not have happed!\n");
		if(DDTR_status>0) printf("There was an error during dedispersion!\n");
		
		
		//------------------> Transfer data to the host
		/*
		if ( (AA_params->enable_acceleration == 1) || (AA_params->enable_periodicity == 1) || (AA_params->analysis_debug ==1) ) {
			for (int k = 0; k < DDTR_plan->ndms[dm_range]; k++) {
				//save_data_offset(d_DDTR_output, (int) (k*DDTR_plan->t_processed[dm_range][t]), h_output_buffer[dm_range][k], ((int) processed_samples)/DDTR_plan->inBin[dm_range], sizeof(float)*((int) DDTR_plan->t_processed[dm_range][t]));
			}
		}
		*/
		// Copy to host would be handled like candidate list in SPS
		//-----------------------------------------------

		if (AA_params->output_dedispersion_data == 1) {
			//for (int k = 0; k < ndms[dm_range]; k++)
			//	write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, output_buffer[dm_range][k], gpu_outputsize, dm_low, dm_high);
		}
			
		checkCudaErrors(cudaGetLastError());
			
		if (AA_params->enable_analysis == 1) {
			int dm_range = DDTR.getCurrentRange();
			int time_chunk = DDTR.getCurrentTimeChunk();
			
			// Set maximum boxcar widths performed between decimations (This should be set at the beginning, so we do not need to set it before every run of SPDT) 
			SPS_params->add_BC_width(PD_MAXTAPS);
			SPS_params->add_BC_width(16);
			SPS_params->add_BC_width(16);
			SPS_params->add_BC_width(16);
			SPS_params->add_BC_width(8);
			SPS_params->print();
			
			// -----------------------------------------------------------------------
			// ------------------------------ New way --------------------------------
			// -----------------------------------------------------------------------
			// TODO:
			// 1) make it possible to use persistent memory allocation for multiple searches, both on the host as well as on the device.
			SPS_Search_AA SPS_search;
			
			SPS_search.setParameters(SPS_params);
			SPS_search.setMSDParameters(MSD_params);
			
			// TODO: DDRT should either return or have a getter for DDTR_OutputData object which would contain all information neccessary for the SPDT to work.
			printf("dm_range: %d; time chunk: %d;\n", dm_range, time_chunk);
			printf("tstart_local: %f; DDTR->tsamp: %f; DDTR->dm_step: %f; DDTR->dm_low: %f; DDTR->dm_high: %f; DDTR->inBin: %d; DDTR->t_processed: %f; DDTR->ndms: %f;\n", tstart_local, DDTR_plan->tsamp, DDTR_plan->dm_step[dm_range], DDTR_plan->dm_low[dm_range], DDTR_plan->dm_high[dm_range], DDTR_plan->inBin[dm_range], DDTR_plan->t_processed[dm_range][time_chunk], DDTR_plan->ndms[dm_range]);
			
			SPS_search.SPS_data.set(tstart_local, DDTR_plan->tsamp, DDTR_plan->dm_step[dm_range], DDTR_plan->dm_low[dm_range], DDTR_plan->dm_high[dm_range], DDTR_plan->inBin[dm_range], DDTR_plan->t_processed[dm_range][time_chunk], DDTR_plan->ndms[dm_range]);
			
			SPS_search.setInputData(d_DDTR_output);
			
			SPS_search.search();
			
			//SPS_search.export_SPSData();
			// or
			//candidatelist->clp.push_back(SPS_search.exportToSubList());
			
			//SPS_CandidatesList* SPS_search.getCandidates
			//vector<SPS_CandidateList*> SPS_Result;
			// rezervovat dopredu 2000;
			//SPS_Results->pushback(SPS_search.getCandidates());
			
			SPS_search.clear();
			// -----------------------------------------------------------------------<
			
			SPS_params->clear_BC_widths();
		}
		printf("INC:\t%zu\n\n\n", processed_samples);
	}
	processed_samples = DDTR.getProcessedSamples() + DDTR_plan->t_processed[0][DDTR_plan->num_tchunks-1];
	tstart_local = ( DDTR_plan->tsamp*processed_samples );

	timer.Stop();
	float time = timer.Elapsed() / 1000;

	printf("\n\n === OVERALL DEDISPERSION THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

	printf("\n(Performed Brute-Force Dedispersion: %g (GPU estimate)",  time);
	printf("\nAmount of telescope time processed: %f", tstart_local);
	processed_samples = DDTR.getProcessedSamples();
	printf("\nNumber of samples processed: %ld", processed_samples);
	printf("\nReal-time speedup factor: %lf", ( tstart_local ) / time);

	cudaFree(d_DDTR_input);
	cudaFree(d_DDTR_output);
	
	#ifdef EXPORT_DD_DATA
		size_t DMs_per_file;
		int *ranges_to_export;
		ranges_to_export = new int[range];
		for(int f=0; f<range; f++) ranges_to_export[f]=1;
		printf("\n\n");
		printf("Exporting dedispersion data...\n");
		DMs_per_file = Calculate_sd_per_file_from_file_size(1000, processed_samples, 1);
		printf("  DM per file: %d;\n", DMs_per_file);
		Export_DD_data(range, output_buffer, processed_samples, ndms, inBin, dm_low, dm_high, dm_step, "DD_data", ranges_to_export, DMs_per_file);
		delete[] ranges_to_export;
	#endif

	double time_processed = ( tstart_local ) / DDTR_plan->tsamp;
	double dm_t_processed = time_processed * DDTR_plan->total_ndms;
	double all_processed = dm_t_processed * DDTR_plan->nchans;
	printf("\nGops based on %.2lf ops per channel per tsamp: %f", NOPS, ( ( NOPS * all_processed ) / ( time ) ) / 1000000000.0);
	int num_reg = SNUMREG;
	float num_threads = DDTR_plan->total_ndms * ( DDTR_plan->t_processed[0][0] ) / ( num_reg );
	float data_size_loaded = ( num_threads * DDTR_plan->nchans * sizeof(ushort) ) / 1000000000;
	float time_in_sec = time;
	float bandwidth = data_size_loaded / time_in_sec;
	printf("\nDevice global memory bandwidth in GB/s: %f", bandwidth);
	printf("\nDevice shared memory bandwidth in GB/s: %f", bandwidth * ( num_reg ));
	float size_gb = ( DDTR_plan->nchans * ( DDTR_plan->t_processed[0][0] ) * sizeof(float) * 8 ) / 1000000000.0;
	printf("\nTelescope data throughput in Gb/s: %f", size_gb / time_in_sec);

	if (AA_params->enable_periodicity == 1) {
		//
		GpuTimer timer;
		timer.Start();
		//
		//GPU_periodicity(range, nsamp, max_ndms, processed_samples, periodicity_sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original, periodicity_nHarmonics, candidate_algorithm, enable_outlier_rejection, OR_sigma_multiplier);
		//
		timer.Stop();
		float time = timer.Elapsed()/1000;
		printf("\n\n === OVERALL PERIODICITY THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Peroidicity Location: %f (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", processed_samples);
		printf("\nReal-time speedup factor: %f", ( tstart_local ) / ( time ));
	}

	if (AA_params->enable_acceleration == 1) {
		// Input needed for fdas is output_buffer which is DDPlan
		// Assumption: gpu memory is free and available
		//
		GpuTimer timer;
		timer.Start();
		//acceleration_fdas(range, nsamp, max_ndms, processed_samples, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original,  FDAS_params, MSD_params);
		//
		timer.Stop();
		float time = timer.Elapsed()/1000;
		printf("\n\n === OVERALL TDAS THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Acceleration Location: %lf (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", processed_samples);
		printf("\nReal-time speedup factor: %lf", ( tstart_local ) / ( time ));
	}
	
}

// TODO: Rename some of the SPS -> SPDT