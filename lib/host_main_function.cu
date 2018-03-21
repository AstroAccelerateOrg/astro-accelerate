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
#include "headers/device_MSD_Configuration.h"
#include "headers/device_MSD.h"
#include "headers/device_MSD_plane_profile.h"

#include "headers/device_SPS_inplace_kernel.h" //Added by KA
#include "headers/device_SPS_inplace.h" //Added by KA
#include "headers/device_SNR_limited.h" //Added by KA
#include "headers/device_SPS_long.h" //Added by KA
#include "headers/device_threshold.h" //Added by KA
#include "headers/device_single_FIR.h" //Added by KA
#include "headers/device_analysis.h" //Added by KA
#include "headers/device_periods.h" //Added by KA
#include "headers/device_peak_find.h" //Added by KA
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
#include "headers/host_MSD_stratagy.h"
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
	int i,
	int t,
	int dm_range,
	int range,
	int enable_debug,
	int enable_analysis,
	int enable_acceleration,
	int enable_output_ffdot_plan,
	int enable_output_fdas_list,
	int enable_periodicity,
	int output_dmt,
	int enable_zero_dm,
	int enable_zero_dm_with_outliers,
	int enable_rfi,
	int enable_sps_baselinenoise,
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
	clock_t start_time,
	int candidate_algorithm,
	int nb_selected_dm,
	float *selected_dm_low,
	float *selected_dm_high,
	int analysis_debug,
	int failsafe,
	float periodicity_sigma_cutoff,
	int periodicity_nHarmonics
	)
{

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
//	printf("\n\n GPU_memory:\t %zu gpu: %zu",gpu_memory/1024/1024,gpu_inputsize);

	checkCudaErrors(cudaGetLastError());
	
	// Allocate memory on host and device.
	printf("\nAllocate memory CPU ...\n");
	allocate_memory_cpu_output(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output,&gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);
	if(enable_debug == 1) debug(5, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);
	printf("\ndone memory CPU ...\n");

	checkCudaErrors(cudaGetLastError());
	
	// Allocate memory on host and device.
	printf("\nAllocate memory GPU ...\n");
	allocate_memory_gpu(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output,
                        &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);
	if(enable_debug == 1) debug(5, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);
	printf("\ndone memory GPU ...\n");

	checkCudaErrors(cudaGetLastError());


        unsigned long int MSD_data_info;
	size_t MSD_profile_size_in_bytes;
        float *d_MSD_workarea = NULL;
	float *d_MSD_interpolated = NULL; 	
	ushort *d_MSD_output_taps = NULL;
	int *gmem_peak = NULL;
	int h_MSD_DIT_width;
	float *h_peak_list;
	float *h_MSD_interpolated, *h_MSD_DIT;
        size_t max_peak_size;
        size_t peak_pos;
        max_peak_size = (size_t) ( max_ndms*t_processed[0][0]/2 );
        h_peak_list   = (float*) malloc(max_peak_size*4*sizeof(float));
//	cudaMallocHost((void **) &h_peak_list, sizeof(float)*max_peak_size*4);

	checkCudaErrors(cudaGetLastError());

	printf("\nStratagy MSD ...\n");
        stratagy_MSD(max_ndms,max_boxcar_width_in_sec, tsamp, t_processed[0][0], &MSD_data_info, &MSD_profile_size_in_bytes, &h_MSD_DIT_width);
	checkCudaErrors(cudaGetLastError());
        
//	printf("Test data: %lld \t MSD_profile: %zu", MSD_data_info, MSD_profile_size_in_bytes);

	printf("\nAllocate memory MSD ...\n");
	allocate_memory_MSD(&d_MSD_workarea, &d_MSD_output_taps, &d_MSD_interpolated, &h_MSD_interpolated, &h_MSD_DIT, &gmem_peak, MSD_data_info, h_MSD_DIT_width, t_processed[0][0], MSD_profile_size_in_bytes);
	checkCudaErrors(cudaGetLastError());


	printf("\n\n GPU_memory:\t %zu gpu: %zu\n\n",gpu_memory/1024/1024,gpu_inputsize);


	cudaStream_t streams[NUM_STREAMS];
	for (int s=0; s < NUM_STREAMS;s++)
		cudaStreamCreate(&streams[s]);

	
	// Clip RFI
	if (enable_rfi) {
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
		checkCudaErrors(cudaGetLastError());

	printf("\nDe-dispersing...");
	GpuTimer timer;
	timer.Start();


	tsamp_original = tsamp;
	maxshift_original = maxshift;
	int t_pos;


	//float *out_tmp;
	//out_tmp = (float *) malloc(( t_processed[0][0] + maxshift ) * max_ndms * sizeof(float));
	//memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));

	// need to add remainder for num_tchunks
	for (t = 0; t < num_tchunks/NUM_STREAMS; t++) {
		for (int s = 0; s < NUM_STREAMS; s++){
	
		t_pos = t*NUM_STREAMS+s;
		printf("\nt_processed:\t%d, %d", t_processed[0][t_pos], t_pos);
		
		checkCudaErrors(cudaGetLastError());

			load_data(-1, inBin, &d_input[(unsigned short)s*(t_processed[0][0]+maxshift)*nchans], &input_buffer[(long int) ( inc * nchans )], t_processed[0][t_pos], maxshift, nchans, dmshifts,streams[s]);
//			load_data(-1, inBin, &d_input[1], &input_buffer[(long int) ( inc * nchans )], t_processed[0][t], maxshift, nchans, dmshifts,streams[s]);

	

		checkCudaErrors(cudaGetLastError());
		
		if (enable_zero_dm) {
			zero_dm(d_input, nchans, t_processed[0][t_pos]+maxshift);
		}
		
		checkCudaErrors(cudaGetLastError());
		
		if (enable_zero_dm_with_outliers) {
			zero_dm_outliers(d_input, nchans, t_processed[0][t_pos]+maxshift);
	 	}
		
		checkCudaErrors(cudaGetLastError());
	
		corner_turn(&d_input[(unsigned short)s*(t_processed[0][0]+maxshift)*nchans], &d_output[s*(t_processed[0][0]+maxshift)*nchans], nchans, t_processed[0][t_pos] + maxshift,streams[s]);
		
	
		checkCudaErrors(cudaGetLastError());
		
		//if (enable_rfi) {
 		//	rfi_gpu(d_input, nchans, t_processed[0][t]+maxshift);
		//}
		
		checkCudaErrors(cudaGetLastError());
		
		int oldBin = 1;
		for (dm_range = 0; dm_range < range; dm_range++) {
			printf("\n\n%f\t%f\t%f\t%d", dm_low[dm_range], dm_high[dm_range], dm_step[dm_range], ndms[dm_range]), fflush(stdout);
			printf("\nAmount of telescope time processed: %f", tstart_local);
			maxshift = maxshift_original / inBin[dm_range];

			checkCudaErrors(cudaGetLastError());
			
			load_data(dm_range, inBin, &d_input[(unsigned short)(s*(t_processed[0][0]+maxshift)*nchans)], &input_buffer[(long int) ( inc * nchans )], t_processed[dm_range][t_pos], maxshift, nchans, dmshifts,streams[s]);
			
			checkCudaErrors(cudaGetLastError());
			
			if (inBin[dm_range] > oldBin) {
				printf("\nBin process.....\n");
				bin_gpu(d_input, d_output, nchans, t_processed[dm_range - 1][t_pos] + maxshift * inBin[dm_range], streams[s]);
				( tsamp ) = ( tsamp ) * 2.0f;
			}
			
			checkCudaErrors(cudaGetLastError());
			
	dedisperse(dm_range, t_processed[dm_range][t_pos], inBin, dmshifts, &d_input[(unsigned short)s*(t_processed[0][0]+maxshift)*nchans], &d_output[s*(t_processed[0][0]+maxshift)*nchans], nchans, ( t_processed[dm_range][t_pos] + maxshift ), maxshift, &tsamp, dm_low, dm_high, dm_step, ndms, nbits, streams[s], failsafe);
//	dedisperse(dm_range, t_processed[dm_range][t], inBin, dmshifts, d_input, d_output, nchans, ( t_processed[dm_range][t] + maxshift ), maxshift, &tsamp, dm_low, dm_high, dm_step, ndms, nbits, failsafe);
	
			checkCudaErrors(cudaGetLastError());
//			cudaStreamSynchronize(streams[s]);

//			if ( (enable_acceleration == 1) || (enable_periodicity == 1) || (analysis_debug ==1) ) {
			if ( (enable_acceleration == 1) || (1 == 1) || (analysis_debug ==1) ) {
				// gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
				//save_data(d_output, out_tmp, gpu_outputsize);

				//#pragma omp parallel for
				for (int k = 0; k < ndms[dm_range]; k++) {
					//memcpy(&output_buffer[dm_range][k][inc / inBin[dm_range]], &out_tmp[k * t_processed[dm_range][t]], sizeof(float) * t_processed[dm_range][t]);
//					cudaHostRegister(output_buffer[dm_range][k],sizeof(float)*(t_processed[dm_range][t_pos]),0);
					save_data_offset(&d_output[s*(t_processed[0][0]+maxshift)*nchans], k * t_processed[dm_range][t_pos], output_buffer[dm_range][k], inc / inBin[dm_range], sizeof(float) * t_processed[dm_range][t_pos], streams[s]);

				}
//				save_data(d_output, &output_buffer[dm_range][0][((long int)inc)/inBin[dm_range]], gpu_outputsize);
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
					gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t_pos] ) * sizeof(float);
					out_tmp = (float *) malloc(( t_processed[0][0] + maxshift ) * max_ndms * sizeof(float));
					memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));
					save_data(d_output, out_tmp, gpu_outputsize);
					analysis_CPU(dm_range, tstart_local, t_processed[dm_range][t_pos], (t_processed[dm_range][t_pos]+maxshift), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, out_tmp,dm_low, dm_high, dm_step, tsamp, max_boxcar_width_in_sec);
					free(out_tmp);
				}
				else {
//					float *h_peak_list;
//					size_t max_peak_size;
//					size_t peak_pos;
//					max_peak_size = (size_t) ( ndms[dm_range]*t_processed[dm_range][t_pos]/2 );
//					h_peak_list   = (float*) malloc(max_peak_size*4*sizeof(float));
//					cudaMallocHost((void **)&h_peak_list,max_peak_size*4*sizeof(float));

					peak_pos=0;
//					analysis_GPU(h_peak_list, &peak_pos, max_peak_size, dm_range, tstart_local, t_processed[dm_range][t_pos], inBin[dm_range], outBin[dm_range], &maxshift, max_ndms, ndms, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, &d_output[s*(t_processed[0][0]+maxshift)*nchans], dm_low, dm_high, dm_step, tsamp, streams[s], candidate_algorithm, enable_sps_baselinenoise, d_MSD_workarea, d_MSD_output_taps, d_MSD_interpolated, h_MSD_DIT, h_MSD_interpolated, gmem_peak, MSD_data_info);
//					analysis_GPU(h_peak_list, &peak_pos, max_peak_size, dm_range, tstart_local, t_processed[dm_range][t_pos], inBin[dm_range], outBin[dm_range], &maxshift, max_ndms, ndms, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, &d_output[s*(t_processed[0][0]+maxshift)*nchans], dm_low, dm_high, dm_step, tsamp, streams[s], candidate_algorithm, enable_sps_baselinenoise, d_MSD_workarea, d_MSD_output_taps, d_MSD_interpolated, gmem_peak, MSD_data_info);

//					int old_s;
//					old_s == (s+1)%NUM_STREAMS;
//					if (s == 0) old_s=1;
//					else old_s=0;
//					cudaStreamSynchronize(streams[old_s]);
//					free(h_peak_list);
//					cudaFreeHost(h_peak_list);
				}

				// This is for testing purposes and should be removed or commented out
				//analysis_CPU(dm_range, tstart_local, t_processed[dm_range][t], (t_processed[dm_range][t]+maxshift), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, out_tmp,dm_low, dm_high, dm_step, tsamp);
			}
			oldBin = inBin[dm_range];
		}

		//memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));

		inc = inc + t_processed[0][t_pos];
		printf("\nINC:\t%ld", inc);
		tstart_local = ( tsamp_original * inc );
		tsamp = tsamp_original;
		maxshift = maxshift_original;
	}
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
	cudaFree(d_MSD_workarea);
	cudaFree(d_MSD_output_taps);
	cudaFree(d_MSD_interpolated);
	cudaFreeHost(input_buffer);
//	cudaFreeHost(h_peak_list);

	checkCudaErrors(cudaGetLastError());
	
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
	printf("\nDevice shared memory bandwidth in GB/s: %f", bandwidth * ( num_reg ));
	float size_gb = ( nchans * ( t_processed[0][0] ) * sizeof(float) * 8 ) / 1000000000.0;
	printf("\nTelescope data throughput in Gb/s: %f", size_gb / time_in_sec);

	checkCudaErrors(cudaGetLastError());

	if (enable_periodicity == 1) {
		//
		GpuTimer timer;
		timer.Start();
		//
//		GPU_periodicity(range, nsamp, max_ndms, inc, periodicity_sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original, periodicity_nHarmonics, candidate_algorithm, enable_sps_baselinenoise, sigma_constant, h_MSD_DIT, h_MSD_interpolated, 0);
		GPU_periodicity(range, nsamp, max_ndms, inc, periodicity_sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original, periodicity_nHarmonics, candidate_algorithm, enable_sps_baselinenoise, sigma_constant, 0);
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
		// acceleration(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original);
		acceleration_fdas(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff,
						  output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original, enable_fdas_custom_fft, enable_fdas_inbin, enable_fdas_norm, sigma_constant, enable_output_ffdot_plan, enable_output_fdas_list);
		//
		timer.Stop();
		float time = timer.Elapsed()/1000;
		printf("\n\n === OVERALL TDAS THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Acceleration Location: %lf (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", inc);
		printf("\nReal-time speedup factor: %lf", ( tstart_local ) / ( time ));
	}

cudaStreamDestroy(streams[NUM_STREAMS-1]);
//	cudaFreeHost(h_MSD_DIT);
//	cudaFreeHost(h_MSD_interpolated);
}
