#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "AstroAccelerate/params.h"
#include "AstroAccelerate/host_periods.h"

#define OUTPUT_TO_FILE 1
#define OUTPUT_TO_LIST 0

void analysis_CPU(int i, float tstart, int t_processed, int nsamp, int nchans, int maxshift, int max_ndms, int *ndms, int *outBin, float cutoff, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp) {

	FILE	*fp_out;
	char	filename[200];

	int	k, dm_count, remaining_time, bin_factor, counter, pos;

	float	start_time;

	unsigned long int j;
	unsigned long int vals=(unsigned long int)(t_processed*ndms[i]);
	
	float mean, stddev, stddev_orig;

	float *exchange_ptr;
	float *binned_output = (float *)malloc(max_ndms*(t_processed)*sizeof(float)/2+1);
	float *binned_output_next = (float *)malloc(max_ndms*(t_processed)*sizeof(float)/4);
	float *output_list = (float *)malloc(vals*sizeof(float));
	int max_list_size=vals/4;

	double	total;

	//printf("\n\n%f\t%f\t%f\t%d", dm_low[i], dm_high[i], dm_step[i], ndms[i]), fflush(stdout);

	// Calculate the total number of values
	

	//chunk=(int)(vals/24);

	//start_time = ((input_increment/nchans)*tsamp);
	start_time=tstart;
	remaining_time = (t_processed);

	sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", start_time, dm_low[i], dm_high[i]);
	if ((fp_out=fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}

	double AA_time_start, AA_time;
	//*************************** TIMING START ************************* 
	AA_time_start=omp_get_wtime();
	// Calculate the mean
	total  = 0.0;
	#pragma omp parallel for default(shared) private(j) reduction(+:total)
	for(j = 0; j < vals; j++) {
		total += (double)output_buffer[j];
	}
	mean = (float)(total/(double)vals);  // Mean for data sample

	// Calculate standard deviation
	total = 0;

	#pragma omp parallel for default(shared) private(j) reduction(+:total)
	for(j = 0; j < vals; j++) {
		total += (double)((output_buffer[j] - mean)*(output_buffer[j] - mean));
	}
	stddev_orig = (float)sqrt(total / (double)vals); // Stddev for data sample

	// Print mean and stddev
	//printf("\nBin: %d, Mean: %f, Stddev: %f", 1, mean, stddev_orig), fflush(stdout);

	// Apply threshold
	pos=0;
	for (dm_count = 0; dm_count < ndms[i]; dm_count++) {
		for(k = 0; k < remaining_time; k++) {
			if((output_buffer[remaining_time*dm_count + k]-mean)/(stddev_orig) >= cutoff) {
				//---------------- Modified ---------------
				if(OUTPUT_TO_FILE)
					fprintf(fp_out, "%f, %f, %f, %d, %d\n", ((float)k)*tsamp+start_time, dm_low[i] + ((float)dm_count)*dm_step[i], (output_buffer[remaining_time*dm_count + k]-mean)/stddev_orig, i, 1);
					//                                                      time       |           DM                            |        SNR                                                   |  | Taps 
				if(OUTPUT_TO_LIST){
					if(pos<max_list_size){
						output_list[4*pos]=((float)k)*tsamp+start_time;
						output_list[4*pos+1]=dm_low[i] + ((float)dm_count)*dm_step[i];
						output_list[4*pos+2]=(output_buffer[remaining_time*dm_count + k]-mean)/stddev_orig;
						output_list[4*pos+3]=1;
						pos++;
					}
				}
				//---------------- Modified ---------------	
			}	
		}
	}

	#pragma omp parallel for private(dm_count,k)
	for (dm_count = 0; dm_count < ndms[i]; dm_count++) {
		int shift=(remaining_time/2)*dm_count;
		int shift2=(remaining_time)*dm_count;
		for(k = 0; k < remaining_time/2; k++) {
			int shift3=2*k+shift2;
			binned_output[shift + k] = ((output_buffer[shift3]) + (output_buffer[shift3 + 1]))*0.5f;
		}
	}

	bin_factor = outBin[i];
	remaining_time=remaining_time/2;
	vals=vals/2;
	counter=1;

	while(bin_factor > 1) {
		// Calculate standard deviation
		total = 0;

		stddev = stddev_orig/((float)sqrt(2.0f*powf(2,(counter-1)))); // Stddev for data sample
		// Print mean and stddev
		//printf("\nBin: %d, Mean: %f, Stddev: %f", (int)powf(2,counter), mean, stddev), fflush(stdout);

		// Apply threshold
		for (dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(k = 0; k < remaining_time; k++) {
				if((binned_output[remaining_time*dm_count + k]-mean)/(stddev) >= cutoff && binned_output[remaining_time*dm_count + k]+(mean)>0.0f) {
					//---------------- Modified ---------------
					if(OUTPUT_TO_FILE)
						fprintf(fp_out, "%f, %f, %f, %d, %d\n", ((float)k)*tsamp+start_time, dm_low[i] + ((float)dm_count)*dm_step[i], (output_buffer[remaining_time*dm_count + k]-mean)/stddev_orig, i, 1);
						//                                                      time       |           DM                            |        SNR                                                   |  | Taps 
					if(OUTPUT_TO_LIST){
						if(pos<max_list_size){
							output_list[4*pos]=((float)k)*tsamp+start_time;
							output_list[4*pos+1]=dm_low[i] + ((float)dm_count)*dm_step[i];
							output_list[4*pos+2]=(output_buffer[remaining_time*dm_count + k]-mean)/stddev_orig;
							output_list[4*pos+3]=1;
							pos++;
						}
					}
					//---------------- Modified ---------------
				}
			}
		}

		#pragma omp parallel for private(dm_count,k) 
		for (dm_count = 0; dm_count < ndms[i]; dm_count++) {
			int shift=(remaining_time/2)*dm_count;
			int shift2=(remaining_time)*dm_count;
			for(k = 0; k < remaining_time/2; k++) {
				int shift3=2*k+shift2;
				binned_output_next[shift + k] = ((binned_output[shift3]) + (binned_output[shift3 + 1]))*0.5f;
			}
		}
	
		remaining_time=remaining_time/2;
		bin_factor=bin_factor/2;
		vals=vals/2;
		exchange_ptr=binned_output;
		binned_output=binned_output_next;
		binned_output_next=exchange_ptr;
		counter++;
	}
	
	AA_time = omp_get_wtime() - AA_time_start;
	//*************************** TIMING END ************************* 
	
	printf("\n====> TIME:%f\n\n",(float) AA_time);
	
	free(binned_output);
	free(binned_output_next);
	free(output_list);
	fclose(fp_out);
}
