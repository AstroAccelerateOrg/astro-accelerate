#include "headers_mains.h"
#include "device_bin.h"
#include "device_init.h"
#include "device_dedisperse.h"
#include "device_dedispersion_kernel.h"
#include "device_load_data.h"
#include "device_corner_turn.h"
#include "device_save_data.h"
#include "host_allocate_memory.h"
#include "host_analysis.h"
#include "host_periods.h"
#include "host_debug.h"
#include "host_get_file_data.h"
#include "host_get_recorded_data.h"
#include "host_get_user_input.h"
#include "host_help.h"
#include "host_rfi.h"
#include "host_stratagy.h"
#include "host_write_file.h"
#include "params.h"


int main(int argc, char *argv[])
{
	// Intenal code variables

	// File pointers
	FILE *fp=NULL;

	// Counters and flags
	int i, t, dm_range;
	int range=0;
	int enable_debug=0;	
	int enable_analysis=0;	
	int enable_periodicity=0;	
	int output_dmt=0;	
	int *inBin=NULL;
	int *outBin=NULL;
	int *ndms=NULL;
	int maxshift=0;
	int max_ndms=0;
	int max_samps=0;
	int num_tchunks=0;
	int total_ndms=0;
	int multi_file=1;
	
	float max_dm=0.0f;

	// Memory sizes and pointers
        size_t inputsize=0;
        size_t outputsize=0;
	size_t gpu_inputsize=0;
	size_t gpu_outputsize=0;
	size_t gpu_memory=0;

        float *input_buffer=NULL;
	float ***output_buffer=NULL;

	float *d_input=NULL;
	float *d_output=NULL;

	float *dmshifts=NULL;

	float *user_dm_low=NULL;
	float *user_dm_high=NULL;
	float *user_dm_step=NULL;
	float *dm_low=NULL;
	float *dm_high=NULL;
	float *dm_step=NULL;

	// Telescope parameters
	int nchans=0;
	int nsamp=0;
	int nbits=0;
	int nsamples=0;
	int nifs=0;
	int **t_processed;

	int nboots = -1;
	int ntrial_bins;
	int navdms=1;
	int nsearch=3;
	
	float aggression=2.5;
	float narrow=0.001f;
	float wide = 0.1f;

	int	maxshift_original;
	double	tsamp_original;

	long int inc=0;

	float tstart=0.0f;
	float tstart_local=0.0f;
	float tsamp=0.0f;
	float fch1=0.0f;
	float foff=0.0f;
	
	// Analysis variables
	float power=2.0f;
	float sigma_cutoff=6.0f;

	// Timing parameters
	double start_time = omp_get_wtime();

	// Users desired de-dispersion strategy. Pick up user defined values from the CLI.
	get_user_input(&fp, argc, argv, &multi_file, &enable_debug, &enable_analysis, &enable_periodicity, &output_dmt, &nboots, &ntrial_bins, &navdms, &narrow, &wide, &aggression, &nsearch, &inBin, &outBin, &power, &sigma_cutoff, &range, 
                       &user_dm_low, &user_dm_high, &user_dm_step);
	if(enable_debug == 1) debug(1, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, 
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	// Initialise the GPU.	
	init_gpu(argc, argv, enable_debug, &gpu_memory);
	if(enable_debug == 1) debug(2, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, 
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);
	

	// Reads telescope parameters from the header of the input file and then counts the number of samples in the input data file.
	get_file_data(&fp, &nchans, &nsamples, &nsamp, &nifs, &nbits, &tsamp, &tstart, &fch1, &foff);
	if(enable_debug == 1) debug(3, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);


	// Calculate the desipersion stratagy.
	stratagy(&maxshift, &max_samps, &num_tchunks, &max_ndms, &total_ndms, &max_dm, power, nchans, nsamp, fch1, foff, tsamp, range, user_dm_low, user_dm_high, user_dm_step,
                 &dm_low, &dm_high, &dm_step, &ndms, &dmshifts, inBin, &t_processed, &gpu_memory);
	if(enable_debug == 1) debug(4, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, 
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);
	
	// Allocate memory on host and device.
	allocate_memory(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output, 
                        &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);
	if(enable_debug == 1) debug(5, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, 
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);


	// Store the recorded telescope data contained in the input filterbank file in the allocated memory.
	get_recorded_data(&fp, nsamp, nchans, nbits, &input_buffer, &inputsize);
     	if(enable_debug == 1) debug(7, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, 
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);
	

	// Specify texture
	struct cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeLinear;
	resDesc.res.linear.devPtr = input_buffer;
	resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
	resDesc.res.linear.desc.x = 16;
	resDesc.res.linear.sizeInBytes = gpu_inputsize;

	// Specify texture object parameters
	struct cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.readMode = cudaReadModeElementType;

	// Create texture object
	cudaTextureObject_t tex;
	cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL);


	// Clip RFI

	//rfi(nsamp, nchans, &input_buffer);
/*	FILE	*fp_o;

	if ((fp_o=fopen("rfi_clipped.dat", "w")) == NULL) {
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}

	for(int t=0; t<nsamp; t++) {
		for(int c = 0; c < nchans; c++) {
			fprintf(fp_o, "%f ", input_buffer[c + t*nchans]);
		}
		fprintf(fp_o, "\n");
	}
*/
	printf("\nDe-dispersing...");
	double start_t, end_t;
	start_t = omp_get_wtime();

	tsamp_original=tsamp;
	maxshift_original=maxshift;

	float *tmp;
	tmp = (float *)malloc((t_processed[0][0]+maxshift)*nchans*sizeof(float));

	float *out_tmp;
	out_tmp = (float *)malloc((t_processed[0][0]+maxshift)*max_ndms*sizeof(float));
	memset(out_tmp,0.0f,t_processed[0][0]+maxshift*max_ndms*sizeof(float));

	for(t=0; t < num_tchunks; t++) {
		printf("\nt_processed:\t%d, %d", t_processed[0][t], t);
//		#pragma omp parallel for
		for(i=0; i < (t_processed[0][t]+maxshift)*nchans; i=i+4) {
			tmp[i  ] = input_buffer[(long int)(inc*nchans)+i];
			tmp[i+1] = input_buffer[(long int)(inc*nchans)+i+1];
			tmp[i+2] = input_buffer[(long int)(inc*nchans)+i+2];
			tmp[i+3] = input_buffer[(long int)(inc*nchans)+i+3];
		}
		//rfi((t_processed[0][t]+maxshift), nchans, &tmp);

		load_data(-1, inBin, d_input, tmp, t_processed[0][t], maxshift, nchans, dmshifts);
		corner_turn(d_input, d_output, nchans, t_processed[0][t]+maxshift); 

		for(dm_range=0; dm_range < range; dm_range++) {
			printf("\n\n%f\t%f\t%f\t%d", dm_low[dm_range], dm_high[dm_range], dm_step[dm_range], ndms[dm_range]), fflush(stdout);
			printf("\nAmount of telescope time processed: %f", tstart_local);
			maxshift=maxshift_original/inBin[dm_range];

			cudaDeviceSynchronize();
			load_data(dm_range, inBin, d_input, tmp, t_processed[dm_range][t], maxshift, nchans, dmshifts);

			if(inBin[dm_range] > 1) {
				bin_gpu(d_input, d_output, nchans, t_processed[dm_range-1][t]+maxshift*inBin[dm_range]);
				(tsamp)=(tsamp)*2.0f;
			}

			dedisperse(dm_range, t_processed[dm_range][t], inBin, dmshifts, d_input, tex, d_output, nchans, (t_processed[dm_range][t]+maxshift), maxshift, &tsamp, dm_low, dm_high, dm_step, ndms);

			gpu_outputsize=ndms[dm_range]*(t_processed[dm_range][t])*sizeof(float);
			cudaDeviceSynchronize();

			save_data(d_output, out_tmp, gpu_outputsize);

			for(int k = 0; k < ndms[dm_range]; k++) {
				for(int l = 0; l < t_processed[dm_range][t]; l++) {
					output_buffer[dm_range][k][l+inc/inBin[dm_range]] = out_tmp[k*t_processed[dm_range][t]+l];
				}					
			}
			if(output_dmt == 1) write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, out_tmp, gpu_outputsize, dm_low, dm_high);		
			if(enable_analysis == 1) analysis(dm_range, tstart_local, t_processed[dm_range][t], (t_processed[dm_range][t]+maxshift), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, out_tmp,
                                                          dm_low, dm_high, dm_step, tsamp);
		}

		memset(out_tmp,0.0f,t_processed[0][0]+maxshift*max_ndms*sizeof(float));

		inc=inc+t_processed[0][t];
		printf("\nINC:\t%d", inc);
		tstart_local=(tsamp_original*inc);
		tsamp=tsamp_original;
		maxshift=maxshift_original;
	}
	end_t=omp_get_wtime();

	printf("\n\n === OVERALL DEDISPERSION THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

	float time = (float)(end_t-start_t);
	printf("\nPerformed Brute-Force Dedispersion: %f (GPU estimate)", time);
	printf("\nAmount of telescope time processed: %f", tstart_local);
	printf("\nNumber of samples processed: %d", inc);
	printf("\nReal-time speedup factor: %f", (tstart_local)/(time));

	cudaFree(d_input);
	cudaFree(d_output);


	double time_processed = (tstart_local)/tsamp_original;
	double dm_t_processed = time_processed * total_ndms;
	double all_processed = dm_t_processed * nchans;
	printf("\nGops based on %.2f ops per channel per tsamp: %f",NOPS,((NOPS*all_processed)/(time))/1000000000.0);
	int num_reg         = SNUMREG;
	float num_threads = total_ndms*(t_processed[0][0])/(num_reg);
	float data_size_loaded = (num_threads * nchans * sizeof(float))/1000000000;
	float time_in_sec = time;
	float bandwidth = data_size_loaded / time_in_sec;
	printf("\nDevice global memory bandwidth in GB/s: %f", bandwidth);
	printf("\nDevice shared memory bandwidth in GB/s: %f", bandwidth*(num_reg));
	float size_gb = (nchans*(t_processed[0][0])*sizeof(float)*8)/1000000000.0;
	printf("\nTelescope data throughput in Gb/s: %f", size_gb/time_in_sec);

	start_t=omp_get_wtime();
	if(enable_periodicity == 1) periodicity(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original);
	end_t=omp_get_wtime();
	time = (float)(end_t-start_t);

	printf("\n\n === OVERALL ACCELERATION THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

	printf("\nPerformed Acceleration Location: %f (GPU estimate)", time);
	printf("\nAmount of telescope time processed: %f", tstart_local);
	printf("\nNumber of samples processed: %d", inc);
	printf("\nReal-time speedup factor: %f", (tstart_local)/(time));

	fclose(fp);

	free(input_buffer);
	free(output_buffer);

	return 0;
}
