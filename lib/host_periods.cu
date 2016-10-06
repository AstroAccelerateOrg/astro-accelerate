#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <math.h>
#include "AstroAccelerate/params.h"

<<<<<<< HEAD
void periodicity(int range, int nsamp, int max_ndms, int processed, int nboots, int num_trial_bins, int navdms, float narrow, float wide, int nsearch, float aggression, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp)
{
	int nharms = 32;

	// Example FFT....
	printf("/n[1DCUFFT] is starting...\n");

	FILE	*fp_c;
	// FILE    *fp_dm, *fp_harm;
	//FILE	*fp_c, *fp_dm;
	char	filename[200];
	
	for(int i=0; i<range; i++)
	{
=======
void periodicity(int range, int nsamp, int max_ndms, int processed, int nboots, int num_trial_bins, int navdms, float narrow, float wide, int nsearch, float aggression, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp) {

	// Example FFT....

	printf("\n");

	printf("[1DCUFFT] is starting...\n");

	FILE	*fp_c, *fp_dm, *fp_harm;
	char	filename[200];
	
	int number_of_candidates=10;

	float** h_top_list = (float**)malloc(sizeof(float*)*5);
	for(int a=0; a < 5; a++) {
		h_top_list[a] = (float*)malloc(sizeof(float)*number_of_candidates);
	}
	for(int a=0; a < 5; a++) {
		for(int b=0; b< number_of_candidates; b++) {
			h_top_list[a][b] = 0.0f;
		}
	}

	for(int i=0; i<range; i++) {
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
		int samps = processed/inBin[i];

		// Allocate memory for signal
		cufftReal* d_signal_in;
		cudaMalloc((void**)&d_signal_in, samps*sizeof(cufftReal));

		cufftComplex* d_signal_out;
		cudaMalloc((void**)&d_signal_out, (samps/2 + 1)*sizeof(cufftComplex));

<<<<<<< HEAD
		cufftComplex* h_signal = (cufftComplex*)malloc((samps)*sizeof(cufftComplex));
		float* h_signal_p = (float*)malloc(sizeof(float) * (samps) * ndms[i]);
		float*** h_harm = (float***)malloc(sizeof(float**) * ndms[i]);
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++)
		{
			h_harm[dm_count] = (float**)malloc(sizeof(float*) * nharms);
			for(int harm=0; harm<nharms; harm++)
				h_harm[dm_count][harm] = (float*)malloc(sizeof(float) * (samps/2 + 1));
		}
			
=======
		cufftComplex* h_signal = (cufftComplex*)malloc((samps/2 + 1)*sizeof(cufftComplex));
		float* h_signal_x = (float*)malloc(sizeof(float) * (samps/2 + 1) * ndms[i]);
		float* h_signal_y = (float*)malloc(sizeof(float) * (samps/2 + 1) * ndms[i]);
		float* h_signal_p = (float*)malloc(sizeof(float) * (samps/2 + 1) * ndms[i]);
		float* h_harm = (float*)malloc(sizeof(float) * (samps/2 + 1) * ndms[i]);
		float* h_signal_inter = (float*)malloc(sizeof(float) * 2*(samps/2 + 1) * ndms[i]);
		
		float** h_candidates = (float**)malloc(sizeof(float*)*ndms[i]);
		for(int a=0; a< ndms[i]; a++) {
			h_candidates[a] = (float*)malloc(sizeof(float)*(samps/2 + 1));
		}
		for(int a=0; a< ndms[i]; a++) {
			for(int b=0; b< samps/2 + 1; b++) {
				h_candidates[a][b] = 0.0f;
			}
		}

		sprintf(filename, "fourier-%d.dat", i);
		if ((fp_c=fopen(filename, "w")) == NULL) {
			fprintf(stderr, "Error opening output file!\n");
			exit(0);
		}
		sprintf(filename, "fourier_inter-%d.dat", i);
		if ((fp_dm=fopen(filename, "w")) == NULL) {
			fprintf(stderr, "Error opening output file!\n");
			exit(0);
		}

>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
		// CUFFT plan
		cufftHandle plan;
		cufftPlan1d(&plan, samps, CUFFT_R2C, 1);

<<<<<<< HEAD
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++)
		{
=======
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {

>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
			cudaMemcpy(d_signal_in, output_buffer[i][dm_count], samps*sizeof(float), cudaMemcpyHostToDevice);

			// Transform signal 
			//printf("\nTransforming dm: %f using cufftExecR2C\n", dm);
			cufftExecR2C(plan, (cufftReal *)d_signal_in, (cufftComplex *)d_signal_out);

			// Copy device memory to host
			cudaMemcpy(h_signal, d_signal_out, sizeof(cufftComplex) * (samps/2 + 1) , cudaMemcpyDeviceToHost);

<<<<<<< HEAD
			h_signal_p[dm_count*(samps/2)] = 0.0;
			#pragma omp parallel for
			for(int j=1;j< samps/2;j++)
			{
				h_signal_p[j+dm_count*(samps/2)] = ((h_signal[j].x*h_signal[j].x + h_signal[j].y*h_signal[j].y));
			}
			#pragma omp parallel for
			for(int harm=1; harm<nharms; harm++)
			{
				for(int j=1; j< (samps/2)/harm;j++)
				{
					float sum=0.0f;
					for(int s=1; s<=harm; s++)
						sum += h_signal_p[(int)((float)j*(float)s)+dm_count*(samps/2)];
					h_harm[dm_count][harm-1][j] = sum/harm; 
				}
=======
			h_signal_p[0+dm_count*(samps/2)] = 0.0; 
			#pragma omp parallel for
			for(int j=1;j< samps/2;j++){
			//	h_signal[j].x = h_signal[j].x-h_signal[0].x;
			//	h_signal[j].y = h_signal[j].y-h_signal[0].y;
				h_signal_x[j+dm_count*(samps/2)] = h_signal[j].x;
				h_signal_y[j+dm_count*(samps/2)] = h_signal[j].y;
				h_signal_p[j+dm_count*(samps/2)]= ((h_signal[j].x*h_signal[j].x + h_signal[j].y*h_signal[j].y));
				h_signal_inter[2*j+dm_count*samps]= h_signal_p[j+dm_count*(samps/2)];
				h_signal_inter[2*j+1+dm_count*samps]=0.616850275*((h_signal[j].x-h_signal[j+1].x)*(h_signal[j].x-h_signal[j+1].x) + (h_signal[j].y-h_signal[j+1].y)*(h_signal[j].y-h_signal[j+1].y));
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
			}
		}

		//Destroy CUFFT context
		cufftDestroy(plan);

		// cleanup memory
		free(h_signal);
		cudaFree(d_signal_in);
		cudaFree(d_signal_out);

<<<<<<< HEAD

		for(int harm=1; harm<nharms; harm++)
		{
			double total = 0.0;
			for (int dm_count = 0; dm_count < ndms[i]; dm_count++)
			{
				for(int j=1;j< (samps/2)/(harm);j++)
					total += ((double)(h_harm[dm_count][harm-1][j]));
			}
			double mean = (total/(double)(((samps/2)/(harm))*ndms[i]));  // Mean for data sample
			
			// Calculate standard deviation
			total = 0.0;
			for (int dm_count = 0; dm_count < ndms[i]; dm_count++)
			{
				for(int j=1;j< (samps/2)/(harm); j++)
					total += (double)((h_harm[dm_count][harm-1][j]-(float)mean)*(h_harm[dm_count][harm-1][j]-(float)mean));
			}
			double stddev = sqrt(abs(total) / (double)((samps/2)/(harm+1)*ndms[i])); // Stddev for data sample

			sprintf(filename, "fourier-%d-%d.dat", i, harm);
			if ((fp_c=fopen(filename, "w")) == NULL)
			{
=======
		double mean, stddev;
	
		double	total=0.0;

		// Calculate the mean
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=0;j< (samps/2);j++){
				total += ((double)(h_signal_p[j+dm_count*(samps/2)]));
			}
		}
		mean = (total/(double)((samps/2)*ndms[i]));  // Mean for data sample

		// Calculate standard deviation
		total = 0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=0;j< (samps/2); j++){
				total += (double)((h_signal_p[j+dm_count*(samps/2)]-(float)mean)*(h_signal_p[j+dm_count*(samps/2)]-(float)mean));
			}
		}
		stddev = sqrt(abs(total) / (double)((samps/2)*ndms[i])); // Stddev for data sample
	
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			float dm = dm_low[i]+dm_step[i]*dm_count;
			for(int j=0;j< (samps/2); j++){
				if((float)((h_signal_p[j+dm_count*(samps/2)]-mean)/stddev) > cutoff) {
					fprintf(fp_c, "\n%f\t%f\t%f\t%f", dm, j*((1.0f/tsamp)/(samps)), (float)(((double)h_signal_p[j+dm_count*(samps/2)]-mean)/stddev), (float)atan(h_signal[j].y/h_signal[j].x));

				}
			}
			fprintf(fp_c, "\n");
			for(int j=0;j< samps; j++){
				if((float)((h_signal_inter[j+dm_count*samps]-mean)/stddev) > cutoff) {
					fprintf(fp_dm, "\n%f\t%lf\t%f\t%f", dm, j*((1.0/tsamp)/(2*samps)), (float)(((double)h_signal_inter[j+dm_count*samps]-mean)/stddev), (float)atan(h_signal[j].y/h_signal[j].x));

				}
			}
			fprintf(fp_dm, "\n");
		}
	
		int harm_max=32;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=0;j< (samps/2);j++){
				h_harm[j+dm_count*(samps/2)] = (h_signal_p[j+dm_count*(samps/2)]);
			}
		}
		int harm=1;
		sprintf(filename, "harmonic-%d-%d.dat", i, harm);
		if ((fp_harm=fopen(filename, "w")) == NULL) {
			fprintf(stderr, "Error opening output file!\n");
			exit(0);
		}

		// Calculate the mean
		total=0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=0;j< (samps/2);j++){
				total += ((double)(h_harm[j+dm_count*(samps/2)]));
			}
		}
		mean = (total/(double)((samps/2)*ndms[i]));  // Mean for data sample

		// Calculate standard deviation
		total = 0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for(int j=0;j< (samps/2); j++){
				total += (double)((h_harm[j+dm_count*(samps/2)]-(float)mean)*(h_harm[j+dm_count*(samps/2)]-(float)mean));
			}
		}
		stddev = sqrt(abs(total) / (double)((samps/2)*ndms[i])); // Stddev for data sample
	
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			float dm = dm_low[i]+dm_step[i]*dm_count;
			for(int j=0;j< samps/2; j++){
				float candidate = (float)(((double)h_harm[j+dm_count*(samps/2)]-mean)/stddev);
				if(candidate > cutoff) {
					fprintf(fp_harm, "\n%f\t%f\t%f", dm, j*1*((1.0f/tsamp)/(samps)), candidate);

					for(int c = 0; c < number_of_candidates; c++) {
						if(candidate > h_top_list[4][c]) {
							for(int d = number_of_candidates - 1; d > c; d--) {
								h_top_list[0][d] = h_top_list[0][d-1];
								h_top_list[1][d] = h_top_list[1][d-1];
								h_top_list[2][d] = h_top_list[2][d-1];
								h_top_list[3][d] = h_top_list[3][d-1];
								h_top_list[4][d] = h_top_list[4][d-1];
							}
							h_top_list[0][c] = dm;
							h_top_list[1][c] = j*1*((1.0f/tsamp)/(samps));
							h_top_list[2][c] = harm;
							h_top_list[3][c] = j;
							h_top_list[4][c] = candidate;
							c=number_of_candidates;
						}
					}
				}
				h_candidates[dm_count][j] = candidate;
			}
			fprintf(fp_harm, "\n");
		}
		fclose(fp_harm);

		for(harm = 2; harm <= harm_max; harm=2*harm) {

			sprintf(filename, "harmonic-%d-%d.dat", i, harm);
			if ((fp_harm=fopen(filename, "w")) == NULL) {
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
				fprintf(stderr, "Error opening output file!\n");
				exit(0);
			}

<<<<<<< HEAD
			for (int dm_count = 0; dm_count < ndms[i]; dm_count++)
			{
				float dm = dm_low[i]+dm_step[i]*dm_count;
				for(int j=1;j< (samps/2)/(harm); j++)
				{
					if((float)(((double)h_harm[dm_count][harm-1][j]-mean)/(stddev/sqrt(harm))) > cutoff)
						fprintf(fp_c, "\n%f\t%f\t%f\t%d", dm, j*((1.0f/tsamp)/(samps)), (float)(((double)h_harm[dm_count][harm-1][j]-mean)/(sqrt(harm)*stddev)), harm);
				}
				fprintf(fp_c, "\n");
			}
			fclose(fp_c);
		}
	}
=======
			for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
				for(int j=0;j < (samps/(2*harm))-harm; j++){
					h_harm[j*harm+dm_count*(samps/2)] += (h_signal_p[j+dm_count*(samps/2)]);
					for(int lerp = j+1; lerp < j+harm; lerp++) h_harm[lerp+dm_count*(samps/2)] += (h_signal_p[j+dm_count*(samps/2)] + 
														      (h_signal_p[j+1+dm_count*(samps/2)]-h_signal_p[j+dm_count*(samps/2)])*((lerp-j)/harm));
				}
			}	
				
			// Calculate the mean
			total=0.0;
			for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
				for(int j=0;j< (samps/2);j++){
					total += ((double)(h_harm[j+dm_count*(samps/2)]));
				}
			}
			mean = (total/(double)((samps/2)*ndms[i]));  // Mean for data sample

			// Calculate standard deviation
			total = 0.0;
			for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
				for(int j=0;j< (samps/2); j++){
					total += (double)((h_harm[j+dm_count*(samps/2)]-(float)mean)*(h_harm[j+dm_count*(samps/2)]-(float)mean));
				}
			}
			stddev = sqrt(abs(total) / (double)((samps/2)*ndms[i])); // Stddev for data sample
			
			for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
				float dm = dm_low[i]+dm_step[i]*dm_count;
				for(int j=0;j< samps/2; j++){
					float candidate = (float)(((double)h_harm[j+dm_count*(samps/2)]-mean)/stddev);
					if(candidate > sqrt(harm)*cutoff) {
						fprintf(fp_harm, "\n%f\t%f\t%f", dm, j*harm*((1.0f/tsamp)/(samps)), candidate);
						for(int c = 0; c < number_of_candidates; c++) {
							if(candidate > h_top_list[4][c]) {
								for(int d = number_of_candidates - 1; d > c; d--) {
									h_top_list[0][d] = h_top_list[0][d-1];
									h_top_list[1][d] = h_top_list[1][d-1];
									h_top_list[2][d] = h_top_list[2][d-1];
									h_top_list[3][d] = h_top_list[3][d-1];
									h_top_list[4][d] = h_top_list[4][d-1];
								}
								h_top_list[0][c] = dm;
								h_top_list[1][c] = j*harm*((1.0f/tsamp)/(samps));
								h_top_list[2][c] = harm;
								h_top_list[3][c] = harm;
								h_top_list[4][c] = candidate;
								c=number_of_candidates;
							}
						}
					}
					h_candidates[dm_count][j] = (float)(((double)h_harm[j+dm_count*(samps/2)]-mean)/stddev);
				}
				fprintf(fp_harm, "\n");
			}
			fclose(fp_harm);
		}
	}
	
	for (int c = 0 ; c < ( number_of_candidates - 1 ); c++) {
		printf("\nCandidate: %d, DM: %f, PERIOD: %f, HARMONIC: %f, PxH: %f, SNR: %f", c, h_top_list[0][c], 1.0f/h_top_list[1][c], h_top_list[2][c], (1.0f/h_top_list[1][c] * h_top_list[2][c]), h_top_list[4][c]);
	}
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
}


