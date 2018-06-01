#include "headers/headers_mains.h"

#include <helper_cuda.h>

#include "headers/device_AA_Parameters.h"

#include "headers/device_bin.h"
#include "headers/device_init.h"
#include "headers/device_dedisperse.h"
#include "headers/device_dedispersion_kernel.h"
#include "headers/device_zero_dm.h"
#include "headers/device_zero_dm_outliers.h"
#include "headers/device_rfi.h"

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



#include "headers/device_load_data.h"
#include "headers/device_corner_turn.h"
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
	
	unsigned short *d_DDTR_input;
	float *d_DDTR_output;
	
	// Allocate memory on host and device.
	allocate_memory_gpu(&d_DDTR_input, &d_DDTR_output, DDTR_plan);
	
	if(AA_params->enable_debug == 1) {
		DDTR_InputData DDTR_data;
		debug(5, start_time, &DDTR_data, DDTR_plan, AA_params, MSD_params, SPS_params, PRS_params, FDAS_params);
	}

	checkCudaErrors(cudaGetLastError());
	
	// Clip RFI
	if (AA_params->enable_rfi) {
		printf("\nPerforming new CPU rfi...");
		rfi(DDTR_plan->nsamp, DDTR_plan->nchans, &h_input_buffer);
	}
	
	printf("\nDe-dispersing...");
	GpuTimer timer;
	timer.Start();


	float tsamp_original = DDTR_plan->tsamp;
	float local_tsamp = DDTR_plan->tsamp;
	int maxshift_original = DDTR_plan->maxshift;
	size_t local_maxshift = DDTR_plan->maxshift;
	float tstart_local = 0;
	size_t inc = 0; //inc = number of processed time samples so far.
	

	for (int t = 0; t < DDTR_plan->num_tchunks; t++) {
		printf("\nt_processed: %d; time chunk: %d; maxshift: %d;\n", DDTR_plan->t_processed[0][t], t, local_maxshift);
		
		checkCudaErrors(cudaGetLastError());

		load_data(-1, DDTR_plan->inBin, d_DDTR_input, &h_input_buffer[inc*DDTR_plan->nchans], (int) DDTR_plan->t_processed[0][t], local_maxshift, (int) DDTR_plan->nchans, DDTR_plan->dmshifts);
		// This will copy the data host->device
		// 		to constant memory it stores dmshifts, nchans, nTimesamples+maxshift and nTimesamples
		// => DDTR.initiateTimeChunk();
		checkCudaErrors(cudaGetLastError());
		
		if (AA_params->enable_zero_dm) {
			zero_dm(d_DDTR_input, (int) DDTR_plan->nchans, (int) (DDTR_plan->t_processed[0][t]+local_maxshift));
		}
		
		checkCudaErrors(cudaGetLastError());
		
		if (AA_params->enable_zero_dm_with_outliers) {
			zero_dm_outliers(d_DDTR_input, (int) DDTR_plan->nchans, (int) (DDTR_plan->t_processed[0][t]+local_maxshift));
	 	}
		
		checkCudaErrors(cudaGetLastError());
	
		corner_turn(d_DDTR_input, d_DDTR_output, (int) DDTR_plan->nchans, (int) (DDTR_plan->t_processed[0][t] + local_maxshift));
		
		checkCudaErrors(cudaGetLastError());
		
		if (AA_params->enable_old_rfi) {
			printf("\nPerforming old GPU rfi...");
 			rfi_gpu(d_DDTR_input, (int) DDTR_plan->nchans, (int)(DDTR_plan->t_processed[0][t]+local_maxshift));
		}
		
		checkCudaErrors(cudaGetLastError());
		
		int oldBin = 1;
		for (int dm_range = 0; dm_range < DDTR_plan->nRanges; dm_range++) {
			printf("\n\n%f\t%f\t%f\t%d", DDTR_plan->dm_low[dm_range], DDTR_plan->dm_high[dm_range], DDTR_plan->dm_step[dm_range], DDTR_plan->ndms[dm_range]), fflush(stdout);
			printf("\nAmount of telescope time processed: %f", tstart_local);
			local_maxshift = maxshift_original/DDTR_plan->inBin[dm_range];

			checkCudaErrors(cudaGetLastError());
			
			cudaDeviceSynchronize();
			
			checkCudaErrors(cudaGetLastError());
			
			load_data(dm_range, DDTR_plan->inBin, d_DDTR_input, &h_input_buffer[inc*DDTR_plan->nchans], (int) DDTR_plan->t_processed[dm_range][t], (int) local_maxshift, (int) DDTR_plan->nchans, DDTR_plan->dmshifts);
			// to constant memory it stores nTimesamples+maxshift and nTimesamples
			// => DDTR.initiateRange(dm_range);
			
			checkCudaErrors(cudaGetLastError());
			
			//---------------------------------
			if (DDTR_plan->inBin[dm_range] > oldBin) {
				bin_gpu(d_DDTR_input, d_DDTR_output, (int) DDTR_plan->nchans, (int) (DDTR_plan->t_processed[dm_range - 1][t] + local_maxshift*DDTR_plan->inBin[dm_range]));
				local_tsamp = local_tsamp*2.0f;
			}
			
			checkCudaErrors(cudaGetLastError());
			
			dedisperse(dm_range, (int) DDTR_plan->t_processed[dm_range][t], DDTR_plan->inBin, DDTR_plan->dmshifts, d_DDTR_input, d_DDTR_output, (int) DDTR_plan->nchans, (int) ( DDTR_plan->t_processed[dm_range][t] + local_maxshift ), (int) local_maxshift, &local_tsamp, DDTR_plan->dm_low, DDTR_plan->dm_high, DDTR_plan->dm_step, DDTR_plan->ndms, DDTR_plan->nbits, AA_params->failsafe);
		
			checkCudaErrors(cudaGetLastError());
			// DDTR.dedisperse();
			//---------------------------------------------
			
			//---------------------------------------------
			if ( (AA_params->enable_acceleration == 1) || (AA_params->enable_periodicity == 1) || (AA_params->analysis_debug ==1) ) {
				for (int k = 0; k < DDTR_plan->ndms[dm_range]; k++) {
					save_data_offset(d_DDTR_output, (int) (k*DDTR_plan->t_processed[dm_range][t]), h_output_buffer[dm_range][k], ((int) inc)/DDTR_plan->inBin[dm_range], sizeof(float)*((int) DDTR_plan->t_processed[dm_range][t]));
				}
			}
			// Copy to host would be handled like candidate list in SPS
			//-----------------------------------------------

			if (AA_params->output_dedispersion_data == 1) {
				//for (int k = 0; k < ndms[dm_range]; k++)
				//	write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, output_buffer[dm_range][k], gpu_outputsize, dm_low, dm_high);
			}
			
			checkCudaErrors(cudaGetLastError());
			
			if (AA_params->enable_analysis == 1) {
				
				printf("\n VALUE OF ANALYSIS DEBUG IS %d\n", AA_params->analysis_debug);

				if (AA_params->analysis_debug == 1) {
					/*
					float *out_tmp;
					gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
					out_tmp = (float *) malloc(( t_processed[0][0] + maxshift ) * max_ndms * sizeof(float));
					memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));
					save_data(d_output, out_tmp, gpu_outputsize);
					analysis_CPU(dm_range, tstart_local, t_processed[dm_range][t], (t_processed[dm_range][t]+maxshift), nchans, local_maxshift, max_ndms, ndms, outBin, sigma_cutoff, out_tmp,dm_low, dm_high, dm_step, tsamp, max_boxcar_width_in_sec);
					free(out_tmp);
					*/
				}
				else {
					// Set maximum boxcar widths performed between decimations
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
					
					
					SPS_search.SPS_data.set(tstart_local, tsamp_original, DDTR_plan->dm_step[dm_range], DDTR_plan->dm_low[dm_range], DDTR_plan->dm_high[dm_range], DDTR_plan->inBin[dm_range], DDTR_plan->t_processed[dm_range][t], DDTR_plan->ndms[dm_range]);
					
					SPS_search.setInputData(d_DDTR_output);
					
					SPS_search.search();
					
					SPS_search.export_SPSData();
					// or
					candidatelist->clp.push_back(SPS_search.exportToSubList());
					
					//SPS_CandidatesList* SPS_search.getCandidates
					//vector<SPS_CandidateList*> SPS_Result;
					// rezervovat dopredu 2000;
					//SPS_Results->pushback(SPS_search.getCandidates());
					
					SPS_search.clear();
					// -----------------------------------------------------------------------<
					
					SPS_params->clear_BC_widths();
				}
			}
			oldBin = DDTR_plan->inBin[dm_range];
		}

		inc = inc + DDTR_plan->t_processed[0][t];
		printf("\nINC:\t%ld", inc);
		tstart_local = ( tsamp_original*inc );
		local_tsamp = tsamp_original;
		local_maxshift = maxshift_original;
	}

	timer.Stop();
	float time = timer.Elapsed() / 1000;

	printf("\n\n === OVERALL DEDISPERSION THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

	printf("\n(Performed Brute-Force Dedispersion: %g (GPU estimate)",  time);
	printf("\nAmount of telescope time processed: %f", tstart_local);
	printf("\nNumber of samples processed: %ld", inc);
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
		DMs_per_file = Calculate_sd_per_file_from_file_size(1000, inc, 1);
		printf("  DM per file: %d;\n", DMs_per_file);
		Export_DD_data(range, output_buffer, inc, ndms, inBin, dm_low, dm_high, dm_step, "DD_data", ranges_to_export, DMs_per_file);
		delete[] ranges_to_export;
	#endif

	double time_processed = ( tstart_local ) / tsamp_original;
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
		//GPU_periodicity(range, nsamp, max_ndms, inc, periodicity_sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original, periodicity_nHarmonics, candidate_algorithm, enable_outlier_rejection, OR_sigma_multiplier);
		//
		timer.Stop();
		float time = timer.Elapsed()/1000;
		printf("\n\n === OVERALL PERIODICITY THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Peroidicity Location: %f (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", inc);
		printf("\nReal-time speedup factor: %f", ( tstart_local ) / ( time ));
	}

	if (AA_params->enable_acceleration == 1) {
		// Input needed for fdas is output_buffer which is DDPlan
		// Assumption: gpu memory is free and available
		//
		GpuTimer timer;
		timer.Start();
		//acceleration_fdas(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original,  FDAS_params, MSD_params);
		//
		timer.Stop();
		float time = timer.Elapsed()/1000;
		printf("\n\n === OVERALL TDAS THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Acceleration Location: %lf (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", inc);
		printf("\nReal-time speedup factor: %lf", ( tstart_local ) / ( time ));
	}
	
}
