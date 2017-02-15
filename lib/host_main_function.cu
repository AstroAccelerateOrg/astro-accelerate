#include "AstroAccelerate/headers_mains.h"
#include "AstroAccelerate/device_bin.h"
#include "AstroAccelerate/device_init.h"
#include "AstroAccelerate/device_dedisperse.h"
#include "AstroAccelerate/device_dedispersion_kernel.h"
#include "AstroAccelerate/device_zero_dm.h"

#include "AstroAccelerate/device_SPS_inplace_kernel.h" //Added by KA
#include "AstroAccelerate/device_SPS_inplace.h" //Added by KA
#include "AstroAccelerate/device_MSD_grid.h" //Added by KA
#include "AstroAccelerate/device_MSD_plane.h" //Added by KA
#include "AstroAccelerate/device_MSD_limited.h" //Added by KA
#include "AstroAccelerate/device_SNR_limited.h" //Added by KA
#include "AstroAccelerate/device_threshold.h" //Added by KA
#include "AstroAccelerate/device_single_FIR.h" //Added by KA
#include "AstroAccelerate/device_analysis.h" //Added by KA

#include "AstroAccelerate/device_load_data.h"
#include "AstroAccelerate/device_corner_turn.h"
#include "AstroAccelerate/device_save_data.h"
#include "AstroAccelerate/host_acceleration.h"
#include "AstroAccelerate/host_allocate_memory.h"
#include "AstroAccelerate/host_analysis.h"
#include "AstroAccelerate/host_periods.h"
#include "AstroAccelerate/host_debug.h"
#include "AstroAccelerate/host_get_file_data.h"
#include "AstroAccelerate/host_get_recorded_data.h"
#include "AstroAccelerate/host_get_user_input.h"
#include "AstroAccelerate/host_help.h"
#include "AstroAccelerate/host_rfi.h"
#include "AstroAccelerate/host_stratagy.h"
#include "AstroAccelerate/host_write_file.h"

#include "AstroAccelerate/host_main_function.h"

#include "AstroAccelerate/params.h"


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
	int enable_rfi,
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
	double start_time
	)
{
	// Initialise the GPU.	
	init_gpu(argc, argv, enable_debug, &gpu_memory);
	if(enable_debug == 1) debug(2, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	// Calculate the desipersion stratagy.
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
	double start_t, end_t;
	start_t = omp_get_wtime();

	tsamp_original = tsamp;
	maxshift_original = maxshift;

	float *out_tmp;
	out_tmp = (float *) malloc(( t_processed[0][0] + maxshift ) * max_ndms * sizeof(float));
	memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));
	

	double start_t_per_range, end_t_per_range, *time_for_range;
	time_for_range = (double *) malloc(range*sizeof(time_for_range));
	for (dm_range=0; dm_range < range; dm_range++) time_for_range[dm_range]=0;

	for (t = 0; t < num_tchunks; t++)
	{
		printf("\nt_processed:\t%d, %d", t_processed[0][t], t);
		//rfi((t_processed[0][t]+maxshift), nchans, &tmp);

		load_data(-1, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[0][t], maxshift, nchans, dmshifts);
		if (enable_zero_dm)
			zero_dm(d_input, nchans, t_processed[0][t]+maxshift);
		corner_turn(d_input, d_output, nchans, t_processed[0][t] + maxshift);
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
			start_t_per_range = omp_get_wtime();
			load_data(dm_range, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[dm_range][t], maxshift, nchans, dmshifts);

			if (inBin[dm_range] > oldBin)
			{
				bin_gpu(d_input, d_output, nchans, t_processed[dm_range - 1][t] + maxshift * inBin[dm_range]);
				( tsamp ) = ( tsamp ) * 2.0f;
			}

			dedisperse(dm_range, t_processed[dm_range][t], inBin, dmshifts, d_input, d_output, nchans, ( t_processed[dm_range][t] + maxshift ), maxshift, &tsamp, dm_low, dm_high, dm_step, ndms);

			gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
			//cudaDeviceSynchronize();

			save_data(d_output, out_tmp, gpu_outputsize);
			//	save_data(d_output, &output_buffer[dm_range][0][((long int)inc)/inBin[dm_range]], gpu_outputsize);

#pragma omp parallel for
			for (int k = 0; k < ndms[dm_range]; k++)
			{
				memcpy(&output_buffer[dm_range][k][inc / inBin[dm_range]], &out_tmp[k * t_processed[dm_range][t]], sizeof(float) * t_processed[dm_range][t]);
			}
			cudaDeviceSynchronize();
			end_t_per_range = omp_get_wtime();
			time_for_range[dm_range] += (end_t_per_range-start_t_per_range);
			if (output_dmt == 1)
				write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, out_tmp, gpu_outputsize, dm_low, dm_high);
			if (enable_analysis == 1) {
				analysis_GPU(dm_range, tstart_local, t_processed[dm_range][t], ( t_processed[dm_range][t] + maxshift ), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, d_output, dm_low, dm_high, dm_step, tsamp);
				
				// This is for testing purposes and should be removed or commented out
				//analysis_CPU(dm_range, tstart_local, t_processed[dm_range][t], (t_processed[dm_range][t]+maxshift), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, out_tmp,dm_low, dm_high, dm_step, tsamp);
				
				
			}
			oldBin = inBin[dm_range];
		}

		memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));

		inc = inc + t_processed[0][t];
		printf("\nINC:\t%ld", inc);
		tstart_local = ( tsamp_original * inc );
		tsamp = tsamp_original;
		maxshift = maxshift_original;
	}
	end_t = omp_get_wtime();

	printf("\n\n === OVERALL DEDISPERSION THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

	float time = (float) ( end_t - start_t );
	printf("\nPerformed Brute-Force Dedispersion: %f (GPU estimate)", time);
	printf("\nAmount of telescope time processed: %f", tstart_local);
	printf("\nNumber of samples processed: %ld", inc);
	printf("\nReal-time speedup factor: %f", ( tstart_local ) / ( time ));

	for (dm_range=0; dm_range<range; dm_range++){
		if (FILTER_OUT_RANGES && dm_range!=RANGE_TO_KEEP) continue;
		printf("\n%d SPEEDUP FACTOR (t processed/sec): %.8f time: %.8f\n", dm_range, tstart_local/time_for_range[dm_range], time_for_range[dm_range]);
	}
	free (time_for_range);

	cudaFree(d_input);
	cudaFree(d_output);
	free(out_tmp);
	free(input_buffer);

	double time_processed = ( tstart_local ) / tsamp_original;
	double dm_t_processed = time_processed * total_ndms;
	double all_processed = dm_t_processed * nchans;
	printf("\nGops based on %.2f ops per channel per tsamp: %f", NOPS, ( ( NOPS * all_processed ) / ( time ) ) / 1000000000.0);
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
		start_t = omp_get_wtime();

		periodicity(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original);

		end_t = omp_get_wtime();
		time = (float) ( end_t - start_t );
		printf("\n\n === OVERALL PERIODICITY THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Peroidicity Location: %f (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", inc);
		printf("\nReal-time speedup factor: %f", ( tstart_local ) / ( time ));
	}

	if (enable_acceleration == 1)
	{
		start_t = omp_get_wtime();

		acceleration(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original);

		end_t = omp_get_wtime();
		time = (float) ( end_t - start_t );
		printf("\n\n === OVERALL TDAS THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Acceleration Location: %f (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", inc);
		printf("\nReal-time speedup factor: %f", ( tstart_local ) / ( time ));
	}

	fclose(fp);

	free(output_buffer);

}
