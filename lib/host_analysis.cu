#include <stdio.h>
#include <stdlib.h>
#include "AstroAccelerate/params.h"
#include "AstroAccelerate/host_periods.h"
#include "AstroAccelerate/device_MSD_plane.h"
#include "AstroAccelerate/device_MSD_limited.h"
#include "AstroAccelerate/device_SNR_limited.h"
#include "AstroAccelerate/device_single_pulse_search.h"
#include "AstroAccelerate/device_threshold.h"
#include "AstroAccelerate/device_single_FIR.h"
#include "timer.h"

//---------------------------------------------------------------------------------
//-------> Kahan MSD
void d_kahan_summation(float *signal, int nDMs, int nTimesamples, int offset, float *result, float *error)
{
	double sum;
	double sum_error;
	double a,b;
	
	sum = 0;
	sum_error = 0;
	for(int d=0; d<nDMs; d++)
	{
		for(int s=0; s<(nTimesamples-offset); s++)
		{
			a = signal[d*nTimesamples + s]-sum_error;
			b = sum+a;
			sum_error = (b-sum);
			sum_error = sum_error-a;
			sum = b;
		}
	}
	*result = sum;
	*error = sum_error;
}

void d_kahan_sd(float *signal, int nDMs, int nTimesamples, int offset, double mean, float *result, float *error)
{
	double sum;
	double sum_error;
	double a,b,dtemp;
	
	sum = 0;
	sum_error = 0;
	for(int d=0; d<nDMs; d++)
	{
		for(int s=0; s<(nTimesamples-offset); s++)
		{
			dtemp = (signal[d*nTimesamples + s]-sum_error - mean);
			a = dtemp*dtemp;
			b = sum+a;
			sum_error = (b-sum);
			sum_error = sum_error-a;
			sum = b;
		}
	}
	*result = sum;
	*error = sum_error;
}

void MSD_Kahan(float *h_input, int nDMs, int nTimesamples, int offset, double *mean, double *sd)
{
	float error, signal_mean, signal_sd;
	int nElements = nDMs*(nTimesamples-offset);
	
	d_kahan_summation(h_input, nDMs, nTimesamples, offset, &signal_mean, &error);
	signal_mean = signal_mean/nElements;
	
	d_kahan_sd(h_input, nDMs, nTimesamples, offset, signal_mean, &signal_sd, &error);
	signal_sd = sqrt(signal_sd/nElements);

	*mean = signal_mean;
	*sd = signal_sd;
}

//-------> Kahan MSD
//---------------------------------------------------------------------------------
void find_min_i_max(float *h_temp, int nDMs, int nTimesamples, int offset, float *max, float *min)
{
	float signal_max, signal_min;
	signal_max = h_temp[0]; 
	signal_min = h_temp[0];
	for(int d=0; d<nDMs; d++)
	{
		for(int s=0; s<(nTimesamples-offset); s++)
		{
			if(h_temp[d*nTimesamples + s]>signal_max) signal_max = h_temp[d*nTimesamples + s];
			if(h_temp[d*nTimesamples + s]<signal_min) signal_min = h_temp[d*nTimesamples + s];
		}
	}
	*max = signal_max;
	*min = signal_min;
}

void print_to_file(float *list, int size, float tsamp, float start_time, float dm_low, float dm_step,  char *filename)
{
	FILE	*file_out;
	
	if ((file_out=fopen(filename, "w")) == NULL)
	{
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}
	
	for (int f=0; f<size; f++)
		fprintf(file_out, "%f, %f, %f, %f\n", list[4*f+1]*tsamp+start_time, dm_low + list[4*f]*dm_step, list[4*f+2], list[4*f+3]);
	
	fclose(file_out);
}

void export_file_nDM_nTimesamples(float *data, int nDMs, int nTimesamples, char *filename)
{
	FILE	*file_out;
	char str[200];
		
	sprintf(str,"%s_DM.dat",filename);
	if ((file_out=fopen(str, "w")) == NULL)
	{
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}
	
	printf("export nDMs\n");
	for (int s=0; s<nTimesamples; s++)
	{
		for(int d=0; d<nDMs; d++)
			fprintf(file_out, "%f ", data[d*nTimesamples + s]);
		fprintf(file_out, "\n");
	}
	
	fclose(file_out);
	
	sprintf(str,"%s_Time.dat",filename);
	if ((file_out=fopen(str, "w")) == NULL)
	{
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}
	
	printf("export nTimesamples\n");
	for(int d=0; d<nDMs; d++)
	{
		for (int s=0; s<nTimesamples; s++)
			fprintf(file_out, "%f ", data[d*nTimesamples + s]);
		fprintf(file_out, "\n");
	}
	
	fclose(file_out);	
}


void analysis(int i, float tstart, int t_processed, int nsamp, int nchans, int maxshift, 
			  int max_ndms, int *ndms, int *outBin, float cutoff, float *output_buffer, 
			  float *dm_low, float *dm_high, float *dm_step, float tsamp)
{
	FILE	*fp_out;
	char	filename[200];

	int remaining_time;

	float	start_time;

	unsigned long int vals;
	int nTimesamples=t_processed;
	int nDMs=ndms[i];
	
	float mean, stddev_orig;

	// Calculate the total number of values
	vals = (unsigned long int)(nDMs*nTimesamples);

	start_time=tstart;
	remaining_time = (t_processed);

	sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", start_time, dm_low[i], dm_high[i]);
	if ((fp_out=fopen(filename, "w")) == NULL)
	{
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}

	float *h_temp=(float*)malloc(vals*sizeof(float));
	float *h_output_list;
	
	//---------------------------------------------------------------------------
	//----------> GPU part
	printf("\n GPU analysis part\n\n");
	printf("Dimensions nDMs:%d; nTimesamples:%d;\n",ndms[i],t_processed);
	GpuTimer timer;
	
	float *d_MSD;			cudaMalloc((void**)&d_MSD, 3*sizeof(float));
	float *d_SNR_MSD;		cudaMalloc((void**)&d_SNR_MSD, 3*sizeof(float));
	float *d_FIR_values;	cudaMalloc((void**)&d_FIR_values, vals*sizeof(float));
	float *d_SNR_values;	cudaMalloc((void**)&d_SNR_values, vals*sizeof(float));		cudaMemset((void*)d_SNR_values, 0 , vals*sizeof(float));
	float *d_SNR_taps;		cudaMalloc((void**)&d_SNR_taps, vals*sizeof(float));		cudaMemset((void*)d_SNR_taps, 0 , vals*sizeof(float));
	int *gmem_pos;			cudaMalloc((void**)&gmem_pos, 1*sizeof(int));				cudaMemset((void*)gmem_pos, 0 , sizeof(int));
	float h_MSD[3];

	int h_list_size;
			
	// ----------------------------------------------------------------------------------------------------	
	// ---------> Mean and standard deviation is calculated once, higher taps are ignored
	MSD_limited(output_buffer, d_MSD, nDMs, nTimesamples, 0);
	
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost);
	mean = h_MSD[0]; 
	stddev_orig = h_MSD[1];
	printf("Bin: %d, Mean: %f, Stddev: %f, nPoints: %f \n", 1, mean, stddev_orig,h_MSD[2]);

	timer.Start();
	PD_SEARCH(output_buffer, d_SNR_values, d_SNR_taps, d_MSD, PD_MAXTAPS-1,  nDMs, nTimesamples);
	timer.Stop();
	printf("PD_SEARCH took:%f ms\n", timer.Elapsed());
	
	timer.Start();
	THRESHOLD_ignore(d_SNR_values, d_SNR_taps, output_buffer, gmem_pos, 10.0, PD_MAXTAPS-1, nDMs, nTimesamples, vals/4);
	timer.Stop();
	printf("THR_WARP took:%f ms\n", timer.Elapsed());
	// ---------> Mean and standard deviation is calculated once, higher taps are ignored
	// ----------------------------------------------------------------------------------------------------	
	
	cudaMemcpy(&h_list_size, gmem_pos, sizeof(int), cudaMemcpyDeviceToHost);
	h_output_list = (float*)malloc(h_list_size*4*sizeof(float));
	cudaMemcpy(h_output_list, output_buffer, h_list_size*4*sizeof(float), cudaMemcpyDeviceToHost);
	
	
	for (int count = 0; count < h_list_size; count++)
		fprintf(fp_out, "%f, %f, %f, %f\n", h_output_list[4*count+1]*tsamp+start_time, dm_low[i] + h_output_list[4*count]*dm_step[i], h_output_list[4*count+2], h_output_list[4*count+3]);
	
	cudaFree(d_MSD);
	cudaFree(d_FIR_values);
	cudaFree(d_SNR_values);
	cudaFree(d_SNR_taps);
	cudaFree(gmem_pos);
	//----------> GPU part
	//---------------------------------------------------------------------------
	
	free(h_temp);
	free(h_output_list);
	
	fclose(fp_out);
}

void export_sps()
{

}