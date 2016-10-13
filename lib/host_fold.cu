#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <math.h>
#include "AstroAccelerate/params.h"

void periodicity(int range, int nsamp, int max_ndms, int processed, int nboots, int num_trial_bins, int navdms, float narrow, float wide, int nsearch, float aggression, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp) {

	int nharms = 32;

	// Example FFT....

	printf("\n");

	printf("[1DCUFFT] is starting...\n");

	FILE	*fp_c, *fp_dm, *fp_harm;
	//FILE	*fp_c, *fp_dm;
	char	filename[200];
	
	for(int i=0; i<range; i++) {
		int samps = processed/inBin[i];

		// Allocate memory for signal
		cufftReal* d_signal_in;
		cudaMalloc((void**)&d_signal_in, samps*sizeof(cufftReal));

		cufftComplex* d_signal_out;
		cudaMalloc((void**)&d_signal_out, (samps/2 + 1)*sizeof(cufftComplex));

		cufftComplex* h_signal = (cufftComplex*)malloc((samps)*sizeof(cufftComplex));
		float* h_signal_p = (float*)malloc(sizeof(float) * (samps) * ndms[i]);
		float*** h_harm = (float***)malloc(sizeof(float**) * ndms[i]);
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			h_harm[dm_count] = (float**)malloc(sizeof(float*) * nharms);
			for(int harm=0; harm<nharms; harm++) {
				h_harm[dm_count][harm] = (float*)malloc(sizeof(float) * (samps/2 + 1));
			}
		}
			
		// CUFFT plan
		cufftHandle plan;
		cufftPlan1d(&plan, samps, CUFFT_R2C, 1);

		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {

			cudaMemcpy(d_signal_in, output_buffer[i][dm_count], samps*sizeof(float), cudaMemcpyHostToDevice);

			// Transform signal 
			//printf("\nTransforming dm: %f using cufftExecR2C\n", dm);
			cufftExecR2C(plan, (cufftReal *)d_signal_in, (cufftComplex *)d_signal_out);

			// Copy device memory to host
			cudaMemcpy(h_signal, d_signal_out, sizeof(cufftComplex) * (samps/2 + 1) , cudaMemcpyDeviceToHost);

			h_signal_p[dm_count*(samps/2)]=0.0;
			#pragma omp parallel for
			for(int j=1;j< samps/2;j++){
				h_signal_p[j+dm_count*(samps/2)]= ((h_signal[j].x*h_signal[j].x + h_signal[j].y*h_signal[j].y));
			}
			#pragma omp parallel for
			for(int harm=1; harm<nharms; harm++) {
				for(int j=1; j< (samps/2)/harm;j++){
					float sum=0.0f;
					for(int s=1; s<=harm; s++) sum+= h_signal_p[(int)((float)j*(float)s)+dm_count*(samps/2)];
					h_harm[dm_count][harm-1][j] = sum/harm; 
				}
			}
		}
/*
		double total = 0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=1;j< (samps/2);j++){
				total += (double)(h_signal_p[j+dm_count*(samps/2)]);
			}
		}
		double mean = (total/(double)(((samps/2))*ndms[i]));  // Mean for data sample
			
		// Calculate standard deviation
		total = 0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=1;j< (samps/2); j++){
				total += (double)((h_signal_p[j+dm_count*(samps/2)]-(float)mean)*(h_signal_p[j+dm_count*(samps/2)]-(float)mean));
			}
		}
		double stddev = sqrt(abs(total) / (double)((samps/2)*ndms[i])); // Stddev for data sample
		printf("mean%lf sd%lf\n", mean, stddev);

		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			h_signal_p[dm_count*(samps/2)]=0.0;
			#pragma omp parallel for
			for(int j=1;j< samps/2;j++){
				h_signal_p[j+dm_count*(samps/2)]= ((h_signal_p[j+dm_count*(samps/2)]+3.0*stddev*((float)rand()/(float)RAND_MAX)) - mean)/stddev;
			}
			cudaMemcpy(d_signal_in, &h_signal_p[dm_count*(samps/2)], samps*sizeof(float), cudaMemcpyHostToDevice);

			// Transform signal 
			//printf("\nTransforming dm: %f using cufftExecR2C\n", dm);
			cufftExecR2C(plan, (cufftReal *)d_signal_in, (cufftComplex *)d_signal_out);

			// Copy device memory to host
			cudaMemcpy(h_signal, d_signal_out, sizeof(cufftComplex) * (samps/2 + 1) , cudaMemcpyDeviceToHost);

			h_signal_p[dm_count*(samps/2)]=0.0;
			#pragma omp parallel for
			for(int j=1;j< samps/2;j++){
				h_signal_p[j+dm_count*(samps/2)]= ((h_signal[j].x*h_signal[j].x + h_signal[j].y*h_signal[j].y));
			}
		}
*/
		//Destroy CUFFT context
		cufftDestroy(plan);

		// cleanup memory
		free(h_signal);
		cudaFree(d_signal_in);
		cudaFree(d_signal_out);
/*
		total = 0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=1;j< (samps/2);j++){
				total += (double)(h_signal_p[j+dm_count*(samps/2)]);
			}
		}
		mean = (total/(double)(((samps/2))*ndms[i]));  // Mean for data sample
			
		// Calculate standard deviation
		total = 0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=1;j< (samps/2); j++){
				total += (double)((h_signal_p[j+dm_count*(samps/2)]-(float)mean)*(h_signal_p[j+dm_count*(samps/2)]-(float)mean));
			}
		}
		stddev = sqrt(abs(total) / (double)((samps/2)*ndms[i])); // Stddev for data sample
		printf("mean%lf sd%lf\n", mean, stddev);

		sprintf(filename, "fourier-%d.dat", i);
		if ((fp_c=fopen(filename, "w")) == NULL) {
			fprintf(stderr, "Error opening output file!\n");
			exit(0);
		}

		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			float dm = dm_low[i]+dm_step[i]*dm_count;
			for(int j=1;j< (samps/2); j++){
				if((h_signal_p[j+dm_count*(samps/2)]-mean)/stddev > cutoff) {
					fprintf(fp_c, "\n%f\t%f\t%f", dm, j*((1.0f/tsamp)/(samps)), (float)(((double)h_signal_p[j+dm_count*(samps/2)])));
				}
			}
			fprintf(fp_c, "\n");
		}
		fclose(fp_c);
*/
/*
		total = 0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=0;j< (samps/2);j++){
				total += (double)(h_signal_p[j+dm_count*(samps/2)]);
			}
		}
		mean = (total/(double)(((samps/2))*ndms[i]));  // Mean for data sample
			
		// Calculate standard deviation
		total = 0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=0;j< (samps/2); j++){
				total += (double)((h_signal_p[j+dm_count*(samps/2)]-(float)mean)*(h_signal_p[j+dm_count*(samps/2)]-(float)mean));
			}
		}
		stddev = sqrt(abs(total) / (double)((samps/2)*ndms[i])); // Stddev for data sample
		printf("mean%f sd%f\n", mean, stddev);

		sprintf(filename, "fourier-%d.dat", i);
		if ((fp_c=fopen(filename, "w")) == NULL) {
			fprintf(stderr, "Error opening output file!\n");
			exit(0);
		}

		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			float dm = dm_low[i]+dm_step[i]*dm_count;
			for(int j=0;j< (samps/2); j++){
				if((float)(((double)h_signal_p[j+dm_count*(samps/2)]-mean)/(stddev)) > cutoff) {
					fprintf(fp_c, "\n%f\t%f\t%f", dm, j*((1.0f/tsamp)/(samps)), (float)(((double)h_signal_p[j+dm_count*(samps/2)]-mean)/(stddev)));
				}
			}
			fprintf(fp_c, "\n");
		}
		fclose(fp_c);
*/


		for(int harm=1; harm<nharms; harm++) {
			double total = 0.0;
			for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
				for(int j=1;j< (samps/2)/(harm);j++){
					total += ((double)(h_harm[dm_count][harm-1][j]));
				}
			}
			double mean = (total/(double)(((samps/2)/(harm))*ndms[i]));  // Mean for data sample
			
			// Calculate standard deviation
			total = 0.0;
			for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
				for(int j=1;j< (samps/2)/(harm); j++){
					total += (double)((h_harm[dm_count][harm-1][j]-(float)mean)*(h_harm[dm_count][harm-1][j]-(float)mean));
				}
			}
			double stddev = sqrt(abs(total) / (double)((samps/2)/(harm+1)*ndms[i])); // Stddev for data sample

			sprintf(filename, "fourier-%d-%d.dat", i, harm);
			if ((fp_c=fopen(filename, "w")) == NULL) {
				fprintf(stderr, "Error opening output file!\n");
				exit(0);
			}

			for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
				float dm = dm_low[i]+dm_step[i]*dm_count;
				for(int j=1;j< (samps/2)/(harm); j++){
					if((float)(((double)h_harm[dm_count][harm-1][j]-mean)/(stddev/sqrt(harm))) > cutoff) {
						fprintf(fp_c, "\n%f\t%f\t%f\t%d", dm, j*((1.0f/tsamp)/(samps)), (float)(((double)h_harm[dm_count][harm-1][j]-mean)/(sqrt(harm)*stddev)), harm);
					}
				}
				fprintf(fp_c, "\n");
			}
			fclose(fp_c);
		}


	}
}


