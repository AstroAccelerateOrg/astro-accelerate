#include <stdio.h>
#include <stdlib.h>
#include "AstroAccelerate/params.h"

#include "AstroAccelerate/device_MSD_plane.h"
#include "AstroAccelerate/device_MSD_limited.h"
#include "AstroAccelerate/device_SNR_limited.h"
#include "AstroAccelerate/device_SPS_inplace.h"
#include "AstroAccelerate/device_threshold.h"
#include "AstroAccelerate/device_single_FIR.h"
#include "timer.h"


void analysis_GPU(int i, float tstart, int t_processed, int nsamp, int nchans, int maxshift, int max_ndms, int *ndms, int *outBin, float cutoff, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp){
	FILE *fp_out;
	char filename[200];
	int max_boxcar_width=16;

	float start_time;
	unsigned long int j;
	unsigned long int vals;
	int nTimesamples = t_processed;
	int nDMs = ndms[i];

	double total;

	// Calculate the total number of values
	vals = (unsigned long int) ( nDMs*nTimesamples );

	//start_time = ((input_increment/nchans)*tsamp);
	start_time = tstart;

	sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", start_time, dm_low[i], dm_high[i]);
	//if ((fp_out=fopen(filename, "w")) == NULL) {
	if (( fp_out = fopen(filename, "wb") ) == NULL)	{
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}

	double signal_mean, signal_sd, total_time, partial_time;
	float signal_mean_1, signal_sd_1, signal_mean_16, signal_sd_16, modifier;
	float max, min, threshold;
	int offset, max_iteration;
	float *h_output_list;

	//---------------------------------------------------------------------------
	//----------> GPU part
	printf("\n GPU analysis part\n\n");
	printf("Dimensions nDMs:%d; nTimesamples:%d;\n", ndms[i], t_processed);
	GpuTimer timer;
	
	float *d_MSD;
	if ( cudaSuccess != cudaMalloc((void**) &d_MSD, sizeof(float)*3)) printf("Allocation error!\n");
	
	float *d_list;
	if ( cudaSuccess != cudaMalloc((void**) &d_list, sizeof(float)*vals)) printf("Allocation error!\n");
	
	float *d_decimated;
	if ( cudaSuccess != cudaMalloc((void **) &d_decimated,  sizeof(float)*(vals/2))) printf("Allocation error!\n");
	
	float *d_boxcar_values;
	if ( cudaSuccess != cudaMalloc((void **) &d_boxcar_values,  sizeof(float)*vals)) printf("Allocation error!\n");
	
	float *d_output_SNR;
	if ( cudaSuccess != cudaMalloc((void **) &d_output_SNR, sizeof(float)*2*vals)) printf("Allocation error!\n");
	
	ushort *d_output_taps;
	if ( cudaSuccess != cudaMalloc((void **) &d_output_taps, sizeof(int)*2*vals)) printf("Allocation error!\n");
	
	
	
	int *gmem_pos;
	cudaMalloc((void**) &gmem_pos, 1*sizeof(int));
	cudaMemset((void*) gmem_pos, 0, sizeof(int));
	float h_MSD[3];

	int h_list_size;

	total_time = 0;
	
	//-------------- Normal MSD
	timer.Start();
	MSD_limited(output_buffer, d_MSD, nDMs, nTimesamples, 128); // Those 128 are there because there was a problem with data, I'm not sure if it is still the case.
	timer.Stop();
	partial_time = timer.Elapsed();
	total_time += partial_time;
	printf("MSD limited took:%f ms\n", partial_time);

	timer.Start();
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost);
	signal_mean_1 = h_MSD[0];
	signal_sd_1 = h_MSD[1];
	printf("MSD Bin: %d, Mean: %f, Stddev: %f\n", 1, signal_mean_1, signal_sd_1);

	offset = PD_FIR(output_buffer, d_list, max_boxcar_width, nDMs, nTimesamples);
	MSD_limited(d_list, d_MSD, nDMs, nTimesamples, offset);
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost);
	signal_mean_16 = h_MSD[0];
	signal_sd_16 = h_MSD[1];
	printf("MSD Bin: %d, Mean: %f, Stddev: %f\n", PD_MAXTAPS, signal_mean_16, signal_sd_16);

	h_MSD[0] = signal_mean_1;
	h_MSD[1] = ( signal_sd_16 - signal_sd_1 )/( (float) ( PD_MAXTAPS - 1 ) );
	h_MSD[2] = signal_sd_1;
	modifier = h_MSD[1];
	cudaMemcpy(d_MSD, h_MSD, 3*sizeof(float), cudaMemcpyHostToDevice);
	timer.Stop();
	partial_time = timer.Elapsed();
	total_time += partial_time;
	printf("Linear sd took:%f ms\n", partial_time);
	//-------------- Normal MSD	
	
	//-------------- SPS	
	timer.Start();
	offset=PD_SEARCH_LONG(output_buffer, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, max_boxcar_width, nDMs, nTimesamples, &max_iteration);
	offset = offset*(1<<(max_iteration-1));
	timer.Stop();
	partial_time = timer.Elapsed();
	total_time += partial_time;
	printf("PD_SEARCH took:%f ms\n", partial_time);
	//-------------- SPS

	//-------------- Thresholding
	timer.Start();
	THRESHOLD(d_output_SNR, d_output_taps, d_list, gmem_pos, cutoff, nDMs, nTimesamples, offset, max_iteration, vals/4);
	timer.Stop();
	partial_time = timer.Elapsed();
	total_time += partial_time;
	printf("THR_WARP took:%f ms\n", partial_time);
	//-------------- Thresholding

	printf("\n====> TOTAL TIME:%f\n\n", total_time);

	cudaMemcpy(&h_list_size, gmem_pos, sizeof(int), cudaMemcpyDeviceToHost);
	h_output_list = (float*) malloc(h_list_size*4*sizeof(float));
	cudaMemcpy(h_output_list, d_list, h_list_size*4*sizeof(float), cudaMemcpyDeviceToHost);

	float *b_list_out;
	b_list_out = (float*) malloc(h_list_size*4*sizeof(float));
	
	#pragma omp parallel for
	for (int count = 0; count < h_list_size; count++){
		b_list_out[4*count] = h_output_list[4*count]*dm_step[i] + dm_low[i];
		b_list_out[4*count + 1] = h_output_list[4*count + 1]*tsamp + start_time;
		b_list_out[4*count + 2] = h_output_list[4*count + 2];
		b_list_out[4*count + 3] = h_output_list[4*count + 3];
	}

	fwrite(&b_list_out[0], h_list_size*sizeof(float), 4, fp_out);

	cudaFree(d_MSD);
	cudaFree(d_list);
	cudaFree(d_boxcar_values);
	cudaFree(d_decimated);
	cudaFree(d_output_SNR);
	cudaFree(d_output_taps);
	cudaFree(gmem_pos);
	//----------> GPU part
	//---------------------------------------------------------------------------

	free(b_list_out);
	free(h_output_list);

	fclose(fp_out);
}
