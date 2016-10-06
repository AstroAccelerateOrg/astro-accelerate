#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
<<<<<<< HEAD
#include <omp.h>
#include "AstroAccelerate/params.h"
#include "AstroAccelerate/device_stats.h"
#include "AstroAccelerate/device_stretch.h"
#include "AstroAccelerate/device_set_stretch.h"
#include "AstroAccelerate/device_power.h"
#include "helper_cuda.h"

void acceleration(int range, int nsamp, int max_ndms, int processed, int nboots, int num_trial_bins, int navdms, float narrow, float wide, int nsearch, float aggression, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp)
{
	// Example FFT....
	printf("\n[1DCUFFT] is starting...\n");

	size_t	size;
	int	a;
	float mean, stddev;

	int chunk = omp_get_num_procs();
	
	for(int i=0; i<range; i++)
	{
		cudaStream_t stream_e;
		cudaError_t result_e;
		result_e = cudaStreamCreate(&stream_e);

		cudaEvent_t event_e;
		cudaEventCreate(&event_e);

		cudaStream_t stream_o;
		cudaError_t result_o;
		result_o = cudaStreamCreate(&stream_o);

		cudaEvent_t event_o;
		cudaEventCreate(&event_o);

		int samps = processed/inBin[i];

		printf("\nsamps:\t%d", samps);
		int nearest = (int)floorf(log2f((float)samps));
		printf("\nnearest:\t%d", nearest);
		samps = (int)powf(2.0, nearest);
		printf("\nsamps:\t%d", samps);


		// Allocate memory for signal even
		float* d_signal_in_e;
		size = samps*sizeof(float);
		printf("\nSize of GPU input signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMalloc((void**)&d_signal_in_e, size));

		float* d_signal_transformed_e;
		size = samps*sizeof(float);
		printf("\nSize of GPU stretched signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMalloc((void**)&d_signal_transformed_e, size));

		cufftComplex* d_signal_fft_e;
		size = (samps/2 + 1)*sizeof(cufftComplex);
		printf("\nSize of GPU output signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMalloc((void**)&d_signal_fft_e, size));

		float* d_signal_power_e;
		size = sizeof(float) * (samps/2) * (2*ACCMAX + ACCSTEP)/ACCSTEP;
		printf("\nSize of GPU power signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMalloc((void**)&d_signal_power_e, size));

		float2* h_signal_e;
		size = (samps)*sizeof(float2);
		printf("\nSize of host output signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMallocHost((void**)&h_signal_e, size));

		float* h_signal_transformed_e;
		size = samps*sizeof(float);
		printf("\nSize of GPU stretched signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMallocHost((void**)&h_signal_transformed_e, size));
	
		float* h_signal_power_e;
		size = sizeof(float) * (samps/2) * (2*ACCMAX + ACCSTEP)/ACCSTEP;
		printf("\nSize of total host power signal:\t%u MB", size/1024/1024), fflush(stdout);
		checkCudaErrors(cudaMallocHost((void**)&h_signal_power_e, size));

		// Allocate memory for signal odd
		float* d_signal_in_o;
		size = samps*sizeof(float);
		printf("\nSize of GPU input signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMalloc((void**)&d_signal_in_o, size));

		float* d_signal_transformed_o;
		size = samps*sizeof(float);
		printf("\nSize of GPU stretched signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMalloc((void**)&d_signal_transformed_o, size));

		cufftComplex* d_signal_fft_o;
		size = (samps/2 + 1)*sizeof(cufftComplex);
		printf("\nSize of GPU output signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMalloc((void**)&d_signal_fft_o, size));

		float* d_signal_power_o;
		size = sizeof(float) * (samps/2) * (2*ACCMAX + ACCSTEP)/ACCSTEP;
		printf("\nSize of GPU power signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMalloc((void**)&d_signal_power_o, size));

		float2* h_signal_o;
		size = (samps)*sizeof(float2);
		printf("\nSize of host output signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMallocHost((void**)&h_signal_o, size));

		float* h_signal_transformed_o;
		size = samps*sizeof(float);
		printf("\nSize of GPU stretched signal:\t%u MB", size/1024/1024);
		checkCudaErrors(cudaMallocHost((void**)&h_signal_transformed_o, size));
	
		float* h_signal_power_o;
		size = sizeof(float) * (samps/2) * (2*ACCMAX + ACCSTEP)/ACCSTEP;
		printf("\nSize of total host power signal:\t%u MB", size/1024/1024), fflush(stdout);
		checkCudaErrors(cudaMallocHost((void**)&h_signal_power_o, size));

		// CUFFT plan even
		cufftHandle plan_e;
		cufftPlan1d(&plan_e, samps, CUFFT_R2C, 1);
		cufftSetStream(plan_e, stream_e);
		
		// CUFFT plan odd
		cufftHandle plan_o;
		cufftPlan1d(&plan_o, samps, CUFFT_R2C, 1);
		cufftSetStream(plan_o, stream_o);
	
		int trials = (2*ACCMAX +ACCSTEP)/ACCSTEP;

		// Transfer even memory asynchronously
		checkCudaErrors(cudaMemcpyAsync(d_signal_in_e, output_buffer[i][0],   samps*sizeof(float), cudaMemcpyHostToDevice, stream_e));
		cudaEventRecord(event_e, stream_e);

		// Cacluclate even dm
		for(a = 0; a < trials; a++)
		{
			int acc = -ACCMAX + a*ACCSTEP;
			float mean = 127.959f;
			set_stretch_gpu(event_e, stream_e, samps, mean, d_signal_transformed_e);
			stretch_gpu(event_e, stream_e, acc, samps, tsamp, d_signal_in_e, d_signal_transformed_e);
			cudaStreamWaitEvent(stream_e, event_e, 0);
			checkCudaErrors(cufftExecR2C(plan_e, (float *)d_signal_transformed_e, (cufftComplex *)d_signal_fft_e));
			power_gpu(event_e, stream_e, samps, a, d_signal_fft_e, d_signal_power_e);
		}

		for (int dm_count = 1; dm_count < ndms[i]-1; dm_count+=2)
		{
			// Transfer odd memory asynchronously
			cudaStreamWaitEvent(stream_o, event_o, 0);
			checkCudaErrors(cudaMemcpyAsync(d_signal_in_o, output_buffer[i][dm_count], samps*sizeof(float), cudaMemcpyHostToDevice, stream_o));
			cudaEventRecord(event_o, stream_o);

			for(a = 0; a < trials; a++) {
				int acc = -ACCMAX + a*ACCSTEP;
				float mean=127.959f;
				set_stretch_gpu(event_o, stream_o, samps, mean, d_signal_transformed_o);
				stretch_gpu(event_o, stream_o, acc, samps, tsamp, d_signal_in_o, d_signal_transformed_o);
				checkCudaErrors(cufftExecR2C(plan_o, (float *)d_signal_transformed_o, (cufftComplex *)d_signal_fft_o));
				cudaStreamWaitEvent(stream_o, event_o, 0);
				power_gpu(event_o, stream_o, samps, a, d_signal_fft_o, d_signal_power_o);
			}

			// Threshold even f-fdot plane
			cudaStreamSynchronize(stream_e);
			stats_gpu(event_e, stream_e, samps, &mean, &stddev, h_signal_power_e, d_signal_power_e);

			// Transfer even memory asynchronously
			checkCudaErrors(cudaMemcpyAsync(d_signal_in_e, output_buffer[i][dm_count+1],   samps*sizeof(float), cudaMemcpyHostToDevice, stream_e));
			cudaEventRecord(event_e, stream_e);
			
			// Cacluclate even dm
			for(a = 0; a < trials; a++) {
				int acc = -ACCMAX + a*ACCSTEP;
				float mean=127.959f;
				set_stretch_gpu(event_e, stream_e, samps, mean, d_signal_transformed_e);
				stretch_gpu(event_e, stream_e, acc, samps, tsamp, d_signal_in_e, d_signal_transformed_e);
				cudaStreamWaitEvent(stream_e, event_e, 0);
				checkCudaErrors(cufftExecR2C(plan_e, (float *)d_signal_transformed_e, (cufftComplex *)d_signal_fft_e));
				power_gpu(event_e, stream_e, samps, a, d_signal_fft_e, d_signal_power_e);
			}

			// Threshold odd f-fdot plane
			cudaStreamSynchronize(stream_o);
			stats_gpu(event_o, stream_o, samps, &mean, &stddev, h_signal_power_o, d_signal_power_o);
		}
		
		//Destroy CUFFT context
		cufftDestroy(plan_e);
		cufftDestroy(plan_o);

		//Destroy streams
		result_e = cudaStreamDestroy(stream_e);
		result_o = cudaStreamDestroy(stream_o);

		// cleanup even memory
		cudaFreeHost(h_signal_e);
		cudaFreeHost(h_signal_power_e);
		cudaFree(d_signal_in_e);
		cudaFree(d_signal_fft_e);
		cudaFree(d_signal_power_e);
		cudaFree(d_signal_transformed_e);

		// cleanup odd memory
		cudaFreeHost(h_signal_o);
		cudaFreeHost(h_signal_power_o);
		cudaFree(d_signal_in_o);
		cudaFree(d_signal_fft_o);
		cudaFree(d_signal_power_o);
		cudaFree(d_signal_transformed_o);
	}
}
=======
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


>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
