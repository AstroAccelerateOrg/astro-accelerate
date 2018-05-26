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
	float **h_SPS_candidatelist,
	size_t *nSPScandidates,
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

	// Initialise the GPU.	
	//init_gpu(argc, argv, enable_debug, &gpu_memory);
	//if(enable_debug == 1) debug(2, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	checkCudaErrors(cudaGetLastError());
	
	// Calculate the dedispersion stratagy.
	//stratagy(&maxshift, &max_samps, &num_tchunks, &max_ndms, &total_ndms, &max_dm, power, nchans, nsamp, fch1, foff, tsamp, range, user_dm_low, user_dm_high, user_dm_step, &dm_low, &dm_high, &dm_step, &ndms, &dmshifts, inBin, &t_processed, &gpu_memory, enable_analysis);
	
	//if(enable_debug == 1) debug(4, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	checkCudaErrors(cudaGetLastError());
	
	// Allocate memory on host and device.
	//allocate_memory_cpu_output(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output, &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);
	//if(enable_debug == 1) debug(5, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	checkCudaErrors(cudaGetLastError());
	
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
	
	/*
	FILE *fp_o;
	if ((fp_o=fopen("rfi_clipped.dat", "wb")) == NULL) {
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}
	fwrite(input_buffer, nchans*nsamp*sizeof(unsigned short), 1, fp_o);
	*/

	printf("\nDe-dispersing...");
	GpuTimer timer;
	timer.Start();


	float tsamp_original = DDTR_plan->tsamp;
	int maxshift_original = DDTR_plan->maxshift;
	float local_tsamp = DDTR_plan->tsamp;
	size_t local_maxshift = DDTR_plan->maxshift;
	float tstart_local = 0;
	size_t inc = 0; //inc = number of processed time samples so far.
	

	for (int t = 0; t < DDTR_plan->num_tchunks; t++) {
		printf("\nt_processed: %d; time chunk: %d; maxshift: %d;\n", DDTR_plan->t_processed[0][t], t, local_maxshift);
		
		checkCudaErrors(cudaGetLastError());

		load_data(-1, DDTR_plan->inBin, d_DDTR_input, &h_input_buffer[inc*DDTR_plan->nchans], (int) DDTR_plan->t_processed[0][t], local_maxshift, (int) DDTR_plan->nchans, DDTR_plan->dmshifts);

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
			
			checkCudaErrors(cudaGetLastError());
			
			if (DDTR_plan->inBin[dm_range] > oldBin) {
				bin_gpu(d_DDTR_input, d_DDTR_output, (int) DDTR_plan->nchans, (int) (DDTR_plan->t_processed[dm_range - 1][t] + local_maxshift*DDTR_plan->inBin[dm_range]));
				local_tsamp = local_tsamp*2.0f;
			}
			
			checkCudaErrors(cudaGetLastError());
			
			dedisperse(dm_range, (int) DDTR_plan->t_processed[dm_range][t], DDTR_plan->inBin, DDTR_plan->dmshifts, d_DDTR_input, d_DDTR_output, (int) DDTR_plan->nchans, (int) ( DDTR_plan->t_processed[dm_range][t] + local_maxshift ), (int) local_maxshift, &local_tsamp, DDTR_plan->dm_low, DDTR_plan->dm_high, DDTR_plan->dm_step, DDTR_plan->ndms, DDTR_plan->nbits, AA_params->failsafe);
		
			checkCudaErrors(cudaGetLastError());
			
			if ( (AA_params->enable_acceleration == 1) || (AA_params->enable_periodicity == 1) || (AA_params->analysis_debug ==1) ) {
				// gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
				//save_data(d_output, out_tmp, gpu_outputsize);

				//#pragma omp parallel for
				for (int k = 0; k < DDTR_plan->ndms[dm_range]; k++) {
					//memcpy(&output_buffer[dm_range][k][inc / inBin[dm_range]], &out_tmp[k * t_processed[dm_range][t]], sizeof(float) * t_processed[dm_range][t]);

					save_data_offset(d_DDTR_output, (int) (k*DDTR_plan->t_processed[dm_range][t]), h_output_buffer[dm_range][k], ((int) inc)/DDTR_plan->inBin[dm_range], sizeof(float)*((int) DDTR_plan->t_processed[dm_range][t]));
				}
			//	save_data(d_output, &output_buffer[dm_range][0][((long int)inc)/inBin[dm_range]], gpu_outputsize);
			}

			if (AA_params->output_dedispersion_data == 1) {
				//for (int k = 0; k < ndms[dm_range]; k++)
				//	write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, output_buffer[dm_range][k], gpu_outputsize, dm_low, dm_high);
				//write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, out_tmp, gpu_outputsize, dm_low, dm_high);
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
					// SPS Parameters test
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
					SPS_Search_AA SPS_search;
					
					SPS_search.SPS_data.set(tstart_local, tsamp_original, DDTR_plan->dm_step[dm_range], DDTR_plan->dm_low[dm_range], DDTR_plan->dm_high[dm_range], DDTR_plan->inBin[dm_range], DDTR_plan->t_processed[dm_range][t], DDTR_plan->ndms[dm_range]);
					
					SPS_search.setParameters(SPS_params);
					SPS_search.setMSDParameters(MSD_params);
					SPS_search.setInputData(d_DDTR_output);
					
					SPS_search.search();
					
					SPS_search.export_SPSData();
					
					SPS_search.clear();
					
					// -----------------------------------------------------------------------<
					
					
					// -----------------------------------------------------------------------
					// ------------------------------ Old way --------------------------------
					// -----------------------------------------------------------------------
					/*
					float *h_peak_list = NULL;
					size_t max_nCandidates = 0;
					size_t peak_pos = 0;
					
					if(SPS_params->candidate_algorithm==1){
						max_nCandidates = (size_t) ( (DDTR_plan->ndms[dm_range]*DDTR_plan->t_processed[dm_range][t])/4 ); // output size for threshold
					}
					else {
						max_nCandidates = (size_t) ( (DDTR_plan->ndms[dm_range]*DDTR_plan->t_processed[dm_range][t])/4 ); // output size for peak-find
					}
					h_peak_list = (float*) malloc(max_nCandidates*4*sizeof(float));
					
					if(h_peak_list==NULL) printf("ERROR: not enough memory to allocate candidate list\n");
					
					// Create description of the data
					SPS_DataDescription SPS_data;
					SPS_data.set(tstart_local, tsamp_original, DDTR_plan->dm_step[dm_range], DDTR_plan->dm_low[dm_range], DDTR_plan->dm_high[dm_range], DDTR_plan->inBin[dm_range], DDTR_plan->t_processed[dm_range][t], DDTR_plan->ndms[dm_range]);
					
					SPS_params->verbose=1;
					analysis_GPU(h_peak_list, &peak_pos, max_nCandidates, SPS_data, d_DDTR_output, SPS_params, MSD_params);
					
					// Current AA behaviour is to change coordinates and write out to disk
					//------------------------> Current AA output
					//#pragma omp parallel for
					//for (int count = 0; count < peak_pos; count++){
					//	h_peak_list[4*count]     = h_peak_list[4*count]*SPS_data.dm_step + SPS_data.dm_low;
					//	h_peak_list[4*count + 1] = h_peak_list[4*count + 1]*SPS_data.sampling_time*SPS_data.inBin + SPS_data.time_start;
					//	h_peak_list[4*count + 2] = h_peak_list[4*count + 2];
					//	h_peak_list[4*count + 3] = h_peak_list[4*count + 3]*SPS_data.inBin;
					//}
					
					char filename[200];
					FILE *fp_out;
					
					if(SPS_params->candidate_algorithm==1){
						sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", SPS_data.time_start, SPS_data.dm_low, SPS_data.dm_high);
					}
					else {
						sprintf(filename, "peak_analysed-t_%.2f-dm_%.2f-%.2f.dat", SPS_data.time_start, SPS_data.dm_low, SPS_data.dm_high);
					}
					
					if(peak_pos>0){
						if (( fp_out = fopen(filename, "wb") ) == NULL)	{
							fprintf(stderr, "Error opening output file!\n");
							exit(0);
						}
						fwrite(h_peak_list, peak_pos*sizeof(float), 4, fp_out);
						fclose(fp_out);
					}
					//------------------------> Current AA output
					free(h_peak_list);
					*/
					// -----------------------------------------------------------------------<
					
					
					
					
					SPS_params->clear_BC_widths();
				}

				// This is for testing purposes and should be removed or commented out
				//analysis_CPU(dm_range, tstart_local, t_processed[dm_range][t], (t_processed[dm_range][t]+maxshift), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, out_tmp,dm_low, dm_high, dm_step, tsamp);
			}
			oldBin = DDTR_plan->inBin[dm_range];
		}

		//memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));

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
