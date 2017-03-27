#include "headers/headers_mains.h"
#include "headers/device_bin.h"
#include "headers/device_init.h"
#include "headers/device_dedisperse.h"
#include "headers/device_dedispersion_kernel.h"
#include "headers/device_zero_dm.h"
#include "headers/device_zero_dm_outliers.h"
#include "headers/device_rfi.h"

#include "headers/device_BLN.h" //Added by KA
#include "headers/device_SPS_inplace_kernel.h" //Added by KA
#include "headers/device_SPS_inplace.h" //Added by KA
#include "headers/device_MSD_grid.h" //Added by KA
#include "headers/device_MSD_plane.h" //Added by KA
#include "headers/device_MSD_limited.h" //Added by KA
#include "headers/device_SNR_limited.h" //Added by KA
#include "headers/device_threshold.h" //Added by KA
#include "headers/device_single_FIR.h" //Added by KA
#include "headers/device_analysis.h" //Added by KA

#include "headers/device_peak_find.h" //Added by KA

#include "headers/device_load_data.h"
#include "headers/device_corner_turn.h"
#include "headers/device_save_data.h"
#include "headers/host_acceleration.h"
#include "headers/host_allocate_memory.h"
#include "headers/host_analysis.h"
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

void main_function
	(
	int argc,
	char* argv[],
	// Internal code variables
	// File pointers
	FILE *fp,
	// Counters and flags
	int i,
	int t,
	int dm_range,
	int range,
	int enable_debug,
	int enable_analysis,
	int enable_acceleration,
	int enable_periodicity,
	int output_dmt,
	int enable_zero_dm,
	int enable_zero_dm_with_outliers,
	int enable_rfi,
	int enable_fdas_custom_fft,
	int enable_fdas_inbin,
	int enable_fdas_norm,
	int *inBin,
	int *outBin,
	int *ndms,
	int maxshift,
	int max_ndms,
	int max_samps,
	int num_tchunks,
	int total_ndms,
	int multi_file,
	float max_dm,
	// Memory sizes and pointers
  size_t inputsize,
  size_t outputsize,
	size_t gpu_inputsize,
	size_t gpu_outputsize,
	size_t gpu_memory,
  unsigned short  *input_buffer,
	float ***output_buffer,
	unsigned short  *d_input,
	float *d_output,
	float *dmshifts,
	float *user_dm_low,
	float *user_dm_high,
	float *user_dm_step,
	float *dm_low,
	float *dm_high,
	float *dm_step,
	// Telescope parameters
	int nchans,
	int nsamp,
	int nbits,
	int nsamples,
	int nifs,
	int **t_processed,
	int nboots,
	int ntrial_bins,
	int navdms,
	int nsearch,
	float aggression,
	float narrow,
	float wide,
	int	maxshift_original,
	double	tsamp_original,
	long int inc,
	float tstart,
	float tstart_local,
	float tsamp,
	float fch1,
	float foff,
	// Analysis variables
	float power,
	float sigma_cutoff,
	float sigma_constant,
	float max_boxcar_width_in_sec,
	clock_t start_time
	)
{
	// Initialise the GPU.	
	init_gpu(argc, argv, enable_debug, &gpu_memory);
	if(enable_debug == 1) debug(2, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	// Calculate the dedispersion stratagy.
	stratagy(&maxshift, &max_samps, &num_tchunks, &max_ndms, &total_ndms, &max_dm, power, nchans, nsamp, fch1, foff, tsamp, range, user_dm_low, user_dm_high, user_dm_step,
                 &dm_low, &dm_high, &dm_step, &ndms, &dmshifts, inBin, &t_processed, &gpu_memory);
	if(enable_debug == 1) debug(4, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	// Allocate memory on host and device.
	allocate_memory_cpu_output(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output,
                        &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);
	if(enable_debug == 1) debug(5, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	// Allocate memory on host and device.
	allocate_memory_gpu(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output,
                        &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);
	if(enable_debug == 1) debug(5, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	// Clip RFI

	//rfi(nsamp, nchans, &input_buffer);
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

	GpuTimer range_timer;
	double *time_for_range;
	time_for_range = (double *) malloc(range*sizeof(time_for_range));
	for (dm_range=0; dm_range < range; dm_range++) time_for_range[dm_range]=0;

	for (t = 0; t < num_tchunks; t++)
	{
		printf("\nt_processed:\t%d, %d", t_processed[0][t], t);

		load_data(-1, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[0][t], maxshift, nchans, dmshifts);

		
		if (enable_zero_dm)
		{
			zero_dm(d_input, nchans, t_processed[0][t]+maxshift);
		}
		
		if (enable_zero_dm_with_outliers)
		{
			zero_dm_outliers(d_input, nchans, t_processed[0][t]+maxshift);
	 	}
	
		corner_turn(d_input, d_output, nchans, t_processed[0][t] + maxshift);
		
		if (enable_rfi)
		{
 			rfi_gpu(d_input, nchans, t_processed[0][t]+maxshift);
		}

		int oldBin = 1;
		for (dm_range = 0; dm_range < range; dm_range++)
		{
			// AB -- this is so that later ranges don't break if an early range is skipped
			if (FILTER_OUT_RANGES && dm_range!=RANGE_TO_KEEP) {
				if (inBin[dm_range] > oldBin)
				{
					bin_gpu(d_input, d_output, nchans, t_processed[dm_range - 1][t] + maxshift * inBin[dm_range]);
					( tsamp ) = ( tsamp ) * 2.0f;
				}
				continue;
			}
			// END AB

			printf("\n\n%f\t%f\t%f\t%d", dm_low[dm_range], dm_high[dm_range], dm_step[dm_range], ndms[dm_range]), fflush(stdout);
			printf("\nAmount of telescope time processed: %f", tstart_local);
			maxshift = maxshift_original / inBin[dm_range];

			cudaDeviceSynchronize();
			range_timer.Start();
			load_data(dm_range, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[dm_range][t], maxshift, nchans, dmshifts);

			if (inBin[dm_range] > oldBin)
			{
				bin_gpu(d_input, d_output, nchans, t_processed[dm_range - 1][t] + maxshift * inBin[dm_range]);
				( tsamp ) = ( tsamp ) * 2.0f;
			}

			dedisperse(dm_range, t_processed[dm_range][t], inBin, dmshifts, d_input, d_output, nchans, ( t_processed[dm_range][t] + maxshift ), maxshift, &tsamp, dm_low, dm_high, dm_step, ndms);

			if (enable_acceleration == 1)
			{
				// gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
				//save_data(d_output, out_tmp, gpu_outputsize);

				//#pragma omp parallel for
				for (int k = 0; k < ndms[dm_range]; k++)
				{
					//memcpy(&output_buffer[dm_range][k][inc / inBin[dm_range]], &out_tmp[k * t_processed[dm_range][t]], sizeof(float) * t_processed[dm_range][t]);

					save_data_offset(d_output, k * t_processed[dm_range][t], output_buffer[dm_range][k], inc / inBin[dm_range], sizeof(float) * t_processed[dm_range][t]);
				}
			//	save_data(d_output, &output_buffer[dm_range][0][((long int)inc)/inBin[dm_range]], gpu_outputsize);
			}
			cudaDeviceSynchronize();
			range_timer.Stop();
			time_for_range[dm_range] += range_timer.Elapsed()/1000.0;;
			if (output_dmt == 1)
			{
				//for (int k = 0; k < ndms[dm_range]; k++)
				//	write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, output_buffer[dm_range][k], gpu_outputsize, dm_low, dm_high);
				//write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, out_tmp, gpu_outputsize, dm_low, dm_high);
			}
			if (enable_analysis == 1) {
				// TODO: put the file export back to analysis I leaving it here at the moment since for interface we need to output from the analysis.
				float *h_output_list;
				float *h_peak_list;
				size_t max_list_size, max_peak_size;
				size_t list_pos, peak_pos;
				max_list_size = (size_t) ( ndms[dm_range]*t_processed[dm_range][t]/2 ); // we can store 1/2 of the input plane
				max_peak_size = (size_t) ( ndms[dm_range]*t_processed[dm_range][t]/2 );
				h_output_list = (float*) malloc(max_list_size*4*sizeof(float)); // Allocations
				h_peak_list   = (float*) malloc(max_list_size*4*sizeof(float));
				
				list_pos=0;
				peak_pos=0;
				
				analysis_GPU(h_output_list, &list_pos, max_list_size, h_peak_list, &peak_pos, max_peak_size, dm_range, tstart_local, t_processed[dm_range][t], inBin[dm_range], outBin[dm_range], &maxshift, max_ndms, ndms, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, d_output, dm_low, dm_high, dm_step, tsamp);
				
				
				printf("-------> list_pos:%zu; \n", list_pos);
				#pragma omp parallel for
				for (int count = 0; count < list_pos; count++){
					h_output_list[4*count]     = h_output_list[4*count]*dm_step[dm_range] + dm_low[dm_range];
					h_output_list[4*count + 1] = h_output_list[4*count + 1]*tsamp + tstart_local;
					//h_output_list[4*count + 2] = h_output_list[4*count + 2];
					//h_output_list[4*count + 3] = h_output_list[4*count + 3];
					
				}
				
				#pragma omp parallel for
				for (int count = 0; count < peak_pos; count++){
					h_peak_list[4*count]     = h_peak_list[4*count]*dm_step[dm_range] + dm_low[dm_range];
					h_peak_list[4*count + 1] = h_peak_list[4*count + 1]*tsamp + tstart_local;
					//h_output_list[4*count + 2] = h_output_list[4*count + 2];
					//h_output_list[4*count + 3] = h_output_list[4*count + 3];
				}

				FILE *fp_out;
				char filename[200];
				
				if(list_pos>0){
					sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", tstart_local, dm_low[dm_range], dm_high[dm_range]);
					//if ((fp_out=fopen(filename, "w")) == NULL) {
					if (( fp_out = fopen(filename, "wb") ) == NULL)	{
						fprintf(stderr, "Error opening output file!\n");
						exit(0);
					}
					fwrite(h_output_list, list_pos*sizeof(float), 4, fp_out);
					fclose(fp_out);
				}
				
				if(peak_pos>0){
					sprintf(filename, "peak_analysed-t_%.2f-dm_%.2f-%.2f.dat", tstart_local, dm_low[dm_range], dm_high[dm_range]);
					//if ((fp_out=fopen(filename, "w")) == NULL) {
					if (( fp_out = fopen(filename, "wb") ) == NULL)	{
						fprintf(stderr, "Error opening output file!\n");
						exit(0);
					}
					fwrite(h_peak_list, peak_pos*sizeof(float), 4, fp_out);
					fclose(fp_out);
				}
				
				
				free(h_peak_list);
				free(h_output_list);
				
				
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

	for (dm_range=0; dm_range<range; dm_range++){
		if (FILTER_OUT_RANGES && dm_range!=RANGE_TO_KEEP) continue;
		printf("\n%d SPEEDUP FACTOR (t processed/sec): %.8f time: %.8f\n", dm_range, tstart_local/time_for_range[dm_range], time_for_range[dm_range]);
	}
	free (time_for_range);

	cudaFree(d_input);
	cudaFree(d_output);
	//free(out_tmp);
	free(input_buffer);

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

	if (enable_periodicity == 1)
	{
		//
		GpuTimer timer;
		timer.Start();
		//
		periodicity(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original);
		//
		timer.Stop();
		float time = timer.Elapsed()/1000;
		printf("\n\n === OVERALL PERIODICITY THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Peroidicity Location: %f (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", inc);
		printf("\nReal-time speedup factor: %f", ( tstart_local ) / ( time ));
	}

	if (enable_acceleration == 1)
	{
		// Input needed for fdas is output_buffer which is DDPlan
		// Assumption: gpu memory is free and available
		//
		GpuTimer timer;
		timer.Start();
		// acceleration(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original);
		acceleration_fdas(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff,
						  output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original, enable_fdas_custom_fft, enable_fdas_inbin, enable_fdas_norm, sigma_constant);
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
