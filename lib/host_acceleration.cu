#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <math.h>
#include "AstroAccelerate/params.h"

void acceleration(int range, int nsamp, int max_ndms, int processed, int nboots, int num_trial_bins, int navdms, float narrow, float wide, int nsearch, float aggression, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp) {

	// Example FFT....

	printf("\n");

	printf("[1DCUFFT] is starting...\n");

	FILE	*fp_c;
	char	filename[200];

	for(int i=0; i<range; i++) {
		int samps = processed/inBin[i];

		// Allocate memory for signal
		cufftReal* d_signal_in;
		cudaMalloc((void**)&d_signal_in, samps*sizeof(cufftReal));

		cufftComplex* d_signal_out;
		cudaMalloc((void**)&d_signal_out, (samps/2 + 1)*sizeof(cufftComplex));

		cufftComplex* h_signal = (cufftComplex*)malloc((samps/2 + 1)*sizeof(cufftComplex));
		float* h_signal_x = (float*)malloc(sizeof(float) * (samps/2 + 1) * ndms[i]);
		float* h_signal_y = (float*)malloc(sizeof(float) * (samps/2 + 1) * ndms[i]);
		float* h_signal_inter_x = (float*)malloc(sizeof(float) * 2*(samps/2 + 1) * ndms[i]);
		float* h_signal_inter_y = (float*)malloc(sizeof(float) * 2*(samps/2 + 1) * ndms[i]);
		
		// CUFFT plan
		cufftHandle plan;
		cufftPlan1d(&plan, samps, CUFFT_R2C, 1);

		sprintf(filename, "acceleration-%d.dat", i);
		if ((fp_c=fopen(filename, "w")) == NULL) {
			fprintf(stderr, "Error opening output file!\n");
			exit(0);
		}
	
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
	
			cudaMemcpy(d_signal_in, output_buffer[i][dm_count], samps*sizeof(float), cudaMemcpyHostToDevice);

			// Transform signal 
			//printf("\nTransforming dm: %f using cufftExecR2C\n", dm);
			cufftExecR2C(plan, (cufftReal *)d_signal_in, (cufftComplex *)d_signal_out);

			// Copy device memory to host
			cudaMemcpy(h_signal, d_signal_out, sizeof(cufftComplex) * (samps/2 + 1) , cudaMemcpyDeviceToHost);

			// Set the DC offset to zero
			//h_signal_p[0+dm_count*(samps/2)] = 0.0; 
			// Store the real and complex parts as floats
			#pragma omp parallel for
			for(int j=1;j< samps/2;j++){
			//	h_signal[j].x = h_signal[j].x-h_signal[0].x;
			//	h_signal[j].y = h_signal[j].y-h_signal[0].y;
				h_signal_x[j+dm_count*(samps/2)] = h_signal[j].x;
				h_signal_y[j+dm_count*(samps/2)] = h_signal[j].y;
				h_signal_inter_x[2*j+dm_count*samps]= h_signal[j+dm_count*(samps/2)].x;
				h_signal_inter_x[2*j+1+dm_count*samps]=0.785398163*((h_signal[j].x-h_signal[j+1].x));
				h_signal_inter_y[2*j+dm_count*samps]= h_signal[j+dm_count*(samps/2)].y;
				h_signal_inter_y[2*j+1+dm_count*samps]=0.785398163*((h_signal[j].y-h_signal[j+1].y));
			}
			int acc_max=0;
			for(int acc=0; acc < acc_max; acc++) {

				// Convolve templates here.....

			}
			for(int j=0;j< samps/2; j++){
				fprintf(fp_c, "\n%d\t%f\t%f", j, h_signal[j].x, h_signal[j].y);
			}
			fprintf(fp_c, "\n");
		}

		//Destroy CUFFT context
		cufftDestroy(plan);

		// cleanup memory
		free(h_signal);
		free(h_signal_x);
		free(h_signal_y);
		free(h_signal_inter_x);
		free(h_signal_inter_y);
		cudaFree(d_signal_in);
		cudaFree(d_signal_out);
	}
}


