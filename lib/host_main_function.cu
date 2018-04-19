#include "headers/headers_mains.h"

#include <helper_cuda.h>

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
#include "headers/device_SPS_Data_Description.h"
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

void main_function
	(
	int argc,
	char* argv[],
	// Internal code variables
	// File pointers
	FILE *fp,
	// Counters and flags
	int i,   //---> AA NOT DOING ANYTHING!
	int t,   //---> AA NOT DOING ANYTHING!
	int dm_range, //---> AA NOT DOING ANYTHING!
	int range, //---> DDTR number of DDTR ranges
	int enable_debug, //---> AA
	int enable_analysis, //---> AA
	int enable_acceleration, //---> AA
	int enable_output_ffdot_plan, //---> FDAS
	int enable_output_fdas_list, //---> FDAS
	int enable_periodicity, //---> AA
	int output_dmt, //---> AA NOT DOING ANYTHING!
	int enable_zero_dm, //---> AA
	int enable_zero_dm_with_outliers, //---> AA
	int enable_rfi, //---> AA
	int enable_old_rfi, //---> AA
	int enable_outlier_rejection, //---> AA WRONG NAME should be MSD but it would be alone
	int enable_fdas_custom_fft, //---> FDAS
	int enable_fdas_inbin, //---> FDAS
	int enable_fdas_norm, //---> FDAS
	int *inBin, //---> DDTR time decimation of data before DDRT for each DDTR range. Defined by user.
	int *outBin, //---> DDTR NOT DOING ANYTHING!
	int *ndms, //---> DDTR number of DM trials in for each DDTR ranges.
	int maxshift, //---> DDTR maxshift calculated by stratagy for given maximum DM search
	int max_ndms, //---> DDTR 
	int max_samps, //---> DDTR NOT DOING ANYTHING!
	int num_tchunks, //---> DDTR number of time chunks into which input data is separated
	int total_ndms, //---> DDTR calculated by stratagy, used for debug!
	int multi_file, //---> AA NOT DOING ANYTHING!
	float max_dm, //---> DDTR calculated by stratagy, used for debug!
	// Memory sizes and pointers
    size_t inputsize, //---> AA I'm not sure about the purpose of this
    size_t outputsize, //---> AA I'm not sure about the purpose of this
	size_t gpu_inputsize, //---> AA I'm not sure about the purpose of this
	size_t gpu_outputsize, //---> AA I'm not sure about the purpose of this
	size_t gpu_memory, //---> AA Free gpu memory (or available memory)
  unsigned short  *input_buffer, //---> AA Input data on the HOST
	float ***output_buffer, //---> AA Output data on the HOST
	unsigned short  *d_input, //---> AA Input data for given time-chunk on the DEVICE
	float *d_output, //---> AA Input data for given time-chunk on the DEVICE
	float *dmshifts, // DDTR internal thing?
	float *user_dm_low,   // DDTR user defined DDTR plan
	float *user_dm_high,  // DDTR user defined DDTR plan
	float *user_dm_step,  // DDTR user defined DDTR plan
	float *dm_low,   // DDTR coordinate transformation
	float *dm_high,  // DDTR coordinate transformation
	float *dm_step,  // DDTR coordinate transformation
	// Telescope parameters
	int nchans, //INPUT/DDTR number of frequency channels
	int nsamp, // INPUT/DDTR number of time-samples
	int nbits, // INPUT/DDTR bit size of input data.
	int nsamples, // INPUT/DDTR should be number of time-samples
	int nifs, // INPUT/DDTR return number of IF channels
	int **t_processed, // DDTR array which contains number of time-samples for each time-chunk from each DM range???
	int nboots, //---> FDAS NOT USED!
	int ntrial_bins, //---> FDAS NOT USED!
	int navdms, //---> FDAS NOT USED!
	int nsearch, //---> FDAS NOT USED!
	float aggression, //---> FDAS NOT USED!
	float narrow, //---> FDAS NOT USED!
	float wide, //---> FDAS NOT USED!
	int	maxshift_original, // NOT NEEDED HERE INTERNAL VARIABLE
	double	tsamp_original, // NOT NEEDED HERE INTERNAL VARIABLE
	size_t inc, // NOT NEEDED HERE INTERNAL VARIABLE total amount of processed timesamples so far.
	float tstart, // INPUT/DDTR does not seem to be used
	float tstart_local, // NOT NEEDED HERE INTERNAL VARIABLE
	float tsamp, // INPUT/DDTR this must not change
	float fch1, // INPUT/DDTR
	float foff, // INPUT/DDTR
	// Analysis variables
	float power, //---> DDTR internal variable of unknown purpose!
	float sigma_cutoff, // SPS/FDAS threshold
	float OR_sigma_multiplier, // MSD outlier rejection rename to OR_sigma_multiplier
	float max_boxcar_width_in_sec, // SPS 
	clock_t start_time, // Time measurements
	int candidate_algorithm, // SPS/FDAS
	int nb_selected_dm, //UNKNOWN
	float *selected_dm_low, //UNKNOWN
	float *selected_dm_high, //UNKNOWN
	int analysis_debug, // AA
	int failsafe, // AA/DDTR
	float periodicity_sigma_cutoff, // PRS
	int periodicity_nHarmonics // PRS
	)
{
	// Configuration of modules
	SPS_Parameters SPS_params(max_boxcar_width_in_sec, sigma_cutoff, candidate_algorithm);
	MSD_Parameters MSD_params(OR_sigma_multiplier, enable_outlier_rejection);
	

	// Initialise the GPU.	
	init_gpu(argc, argv, enable_debug, &gpu_memory);
	if(enable_debug == 1) debug(2, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	checkCudaErrors(cudaGetLastError());
	
	// Calculate the dedispersion stratagy.
	stratagy(&maxshift, &max_samps, &num_tchunks, &max_ndms, &total_ndms, &max_dm, power, nchans, nsamp, fch1, foff, tsamp, range, user_dm_low, user_dm_high, user_dm_step,
                 &dm_low, &dm_high, &dm_step, &ndms, &dmshifts, inBin, &t_processed, &gpu_memory, enable_analysis);
	if(enable_debug == 1) debug(4, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	checkCudaErrors(cudaGetLastError());
	
	// Allocate memory on host and device.
	allocate_memory_cpu_output(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output,
                        &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);
	if(enable_debug == 1) debug(5, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	checkCudaErrors(cudaGetLastError());
	
	// Allocate memory on host and device.
	allocate_memory_gpu(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output,
                        &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);
	if(enable_debug == 1) debug(5, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	checkCudaErrors(cudaGetLastError());
	
	// Clip RFI
	if (enable_rfi) {
		printf("\nPerforming new CPU rfi...");
		rfi(nsamp, nchans, &input_buffer);
	}
	/*
	 FILE	*fp_o;

	 if ((fp_o=fopen("rfi_clipped.dat", "wb")) == NULL) {
	 fprintf(stderr, "Error opening output file!\n");
	 exit(0);
	 }
	 fwrite(input_buffer, nchans*nsamp*sizeof(unsigned short), 1, fp_o);
	 */

	printf("\nDe-dispersing...");
	GpuTimer timer;
	timer.Start();


	tsamp_original = tsamp;
	maxshift_original = maxshift;

	//float *out_tmp;
	//out_tmp = (float *) malloc(( t_processed[0][0] + maxshift ) * max_ndms * sizeof(float));
	//memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));

	for (t = 0; t < num_tchunks; t++) {
		printf("\nt_processed:\t%d, %d", t_processed[0][t], t);
		
		checkCudaErrors(cudaGetLastError());

		load_data(-1, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[0][t], maxshift, nchans, dmshifts);

		checkCudaErrors(cudaGetLastError());
		
		if (enable_zero_dm) {
			zero_dm(d_input, nchans, t_processed[0][t]+maxshift);
		}
		
		checkCudaErrors(cudaGetLastError());
		
		if (enable_zero_dm_with_outliers) {
			zero_dm_outliers(d_input, nchans, t_processed[0][t]+maxshift);
	 	}
		
		checkCudaErrors(cudaGetLastError());
	
		corner_turn(d_input, d_output, nchans, t_processed[0][t] + maxshift);
		
		checkCudaErrors(cudaGetLastError());
		
		if (enable_old_rfi) {
			printf("\nPerforming old GPU rfi...");
 			rfi_gpu(d_input, nchans, t_processed[0][t]+maxshift);
		}
		
		checkCudaErrors(cudaGetLastError());
		
		int oldBin = 1;
		for (dm_range = 0; dm_range < range; dm_range++) {
			printf("\n\n%f\t%f\t%f\t%d", dm_low[dm_range], dm_high[dm_range], dm_step[dm_range], ndms[dm_range]), fflush(stdout);
			printf("\nAmount of telescope time processed: %f", tstart_local);
			maxshift = maxshift_original / inBin[dm_range];

			checkCudaErrors(cudaGetLastError());
			
			cudaDeviceSynchronize();
			
			checkCudaErrors(cudaGetLastError());
			
			load_data(dm_range, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[dm_range][t], maxshift, nchans, dmshifts);
			
			checkCudaErrors(cudaGetLastError());
			
			if (inBin[dm_range] > oldBin) {
				bin_gpu(d_input, d_output, nchans, t_processed[dm_range - 1][t] + maxshift * inBin[dm_range]);
				( tsamp ) = ( tsamp ) * 2.0f;
			}
			
			checkCudaErrors(cudaGetLastError());
			
			dedisperse(dm_range, t_processed[dm_range][t], inBin, dmshifts, d_input, d_output, nchans, ( t_processed[dm_range][t] + maxshift ), maxshift, &tsamp, dm_low, dm_high, dm_step, ndms, nbits, failsafe);
		
			checkCudaErrors(cudaGetLastError());
			
			if ( (enable_acceleration == 1) || (enable_periodicity == 1) || (analysis_debug ==1) ) {
				// gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
				//save_data(d_output, out_tmp, gpu_outputsize);

				//#pragma omp parallel for
				for (int k = 0; k < ndms[dm_range]; k++) {
					//memcpy(&output_buffer[dm_range][k][inc / inBin[dm_range]], &out_tmp[k * t_processed[dm_range][t]], sizeof(float) * t_processed[dm_range][t]);

					save_data_offset(d_output, k * t_processed[dm_range][t], output_buffer[dm_range][k], inc / inBin[dm_range], sizeof(float) * t_processed[dm_range][t]);
				}
			//	save_data(d_output, &output_buffer[dm_range][0][((long int)inc)/inBin[dm_range]], gpu_outputsize);
			}

			if (output_dmt == 1)
			{
				//for (int k = 0; k < ndms[dm_range]; k++)
				//	write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, output_buffer[dm_range][k], gpu_outputsize, dm_low, dm_high);
				//write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, out_tmp, gpu_outputsize, dm_low, dm_high);
			}
			
			checkCudaErrors(cudaGetLastError());
			
			if (enable_analysis == 1) {
				
				printf("\n VALUE OF ANALYSIS DEBUG IS %d\n", analysis_debug);

				if (analysis_debug == 1) {
					float *out_tmp;
					gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
					out_tmp = (float *) malloc(( t_processed[0][0] + maxshift ) * max_ndms * sizeof(float));
					memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));
					save_data(d_output, out_tmp, gpu_outputsize);
					analysis_CPU(dm_range, tstart_local, t_processed[dm_range][t], (t_processed[dm_range][t]+maxshift), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, out_tmp,dm_low, dm_high, dm_step, tsamp, max_boxcar_width_in_sec);
					free(out_tmp);
				}
				else {
					float *h_peak_list;
					size_t max_peak_size;
					size_t peak_pos;
					max_peak_size = (size_t) ( ndms[dm_range]*t_processed[dm_range][t]/2 );
					h_peak_list   = (float*) malloc(max_peak_size*4*sizeof(float));
					
					SPS_Data_Description SPS_data(tstart_local, tsamp_original, dm_step[dm_range], dm_low[dm_range], dm_high[dm_range], inBin[dm_range], t_processed[dm_range][t], ndms[dm_range]);

					peak_pos=0;
					analysis_GPU(h_peak_list, &peak_pos, max_peak_size, SPS_data, d_output, SPS_params, MSD_params);

					free(h_peak_list);
				}

				// This is for testing purposes and should be removed or commented out
				//analysis_CPU(dm_range, tstart_local, t_processed[dm_range][t], (t_processed[dm_range][t]+maxshift), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, out_tmp,dm_low, dm_high, dm_step, tsamp);
			}
			oldBin = inBin[dm_range];
		}

		//memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));

		inc = inc + t_processed[0][t];
		printf("\nINC:\t%ld", inc);
		tstart_local = ( tsamp_original * inc );
		tsamp = tsamp_original;
		maxshift = maxshift_original;
	}

	timer.Stop();
	float time = timer.Elapsed() / 1000;

	printf("\n\n === OVERALL DEDISPERSION THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

	printf("\n(Performed Brute-Force Dedispersion: %g (GPU estimate)",  time);
	printf("\nAmount of telescope time processed: %f", tstart_local);
	printf("\nNumber of samples processed: %ld", inc);
	printf("\nReal-time speedup factor: %lf", ( tstart_local ) / time);

	cudaFree(d_input);
	cudaFree(d_output);
	//free(out_tmp);
	free(input_buffer);
	
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
	double dm_t_processed = time_processed * total_ndms;
	double all_processed = dm_t_processed * nchans;
	printf("\nGops based on %.2lf ops per channel per tsamp: %f", NOPS, ( ( NOPS * all_processed ) / ( time ) ) / 1000000000.0);
	int num_reg = SNUMREG;
	float num_threads = total_ndms * ( t_processed[0][0] ) / ( num_reg );
	float data_size_loaded = ( num_threads * nchans * sizeof(ushort) ) / 1000000000;
	float time_in_sec = time;
	float bandwidth = data_size_loaded / time_in_sec;
	printf("\nDevice global memory bandwidth in GB/s: %f", bandwidth);
	printf("\nDevice shared memory bandwidth in GB/s: %f", bandwidth * ( num_reg ));
	float size_gb = ( nchans * ( t_processed[0][0] ) * sizeof(float) * 8 ) / 1000000000.0;
	printf("\nTelescope data throughput in Gb/s: %f", size_gb / time_in_sec);

	if (enable_periodicity == 1) {
		//
		GpuTimer timer;
		timer.Start();
		//
		GPU_periodicity(range, nsamp, max_ndms, inc, periodicity_sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original, periodicity_nHarmonics, candidate_algorithm, enable_outlier_rejection, OR_sigma_multiplier);
		//
		timer.Stop();
		float time = timer.Elapsed()/1000;
		printf("\n\n === OVERALL PERIODICITY THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Peroidicity Location: %f (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", inc);
		printf("\nReal-time speedup factor: %f", ( tstart_local ) / ( time ));
	}

	if (enable_acceleration == 1) {
		// Input needed for fdas is output_buffer which is DDPlan
		// Assumption: gpu memory is free and available
		//
		GpuTimer timer;
		timer.Start();
		acceleration_fdas(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff,
						  output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original, enable_fdas_custom_fft, enable_fdas_inbin, enable_fdas_norm, OR_sigma_multiplier, enable_output_ffdot_plan, enable_output_fdas_list);
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
