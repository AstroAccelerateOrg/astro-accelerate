#ifndef ZERODM_OUTLIERS_KERNEL_H_
#define ZERODM_OUTLIERS_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

#define MEAN	127.5f
#define CUT 	2.0f
#define R_CUT	4.0f
#define ITER 	20
#define ACC	0.000001f

//{{{ Dristribution based RFI removal - needs cleaning and optimizing // WA 07/08/18
__global__ void zero_dm_outliers_kernel(unsigned short *d_input, int nchans, int nsamp)
{

	int t  = blockIdx.x * blockDim.x + threadIdx.x;

	int count = 0;
	int iters = 0;

	float stdev = 1000000.0f;
	float mean = MEAN;
	float mean_last = 0.0f;
	float sum = 0.0f;
	float sum_squares = 0.0f;
	float cutoff = (CUT * stdev); 

	while(abs(mean - mean_last) > ACC) {
		sum = 0.0f;
		sum_squares = 0.0f;
		count = 0;

		for(int c = 0; c < nchans; c++) {
			float data=(float)d_input[t*nchans + c];
			if(data < (mean + cutoff) && data > (mean - cutoff) ) {
				sum+=data;
				sum_squares+=(data*data);
				count++;
			}
		}
		mean_last = mean;
		mean = (sum/(float)count);
		sum_squares = ((sum_squares / count) - (mean * mean));
	    	stdev = sqrt(sum_squares);
		cutoff = (CUT * stdev); 

		iters++;
		if(iters > ITER) break;
	}

	for(int c = 0; c < nchans; c++) {
		if((d_input[t*nchans + c]-mean < R_CUT*stdev) && (d_input[t*nchans + c]-mean > - R_CUT*stdev)) {
			d_input[t*nchans + c] = (unsigned short)((float)d_input[t*nchans + c]-(float)mean+MEAN);
			//d_input[t*nchans + c] = (unsigned short)((float)d_input[t*nchans + c]-(float)mean_of_mean+MEAN);
		} else {
			//d_input[t*nchans + c] = mean;
			//d_input[t*nchans + c] = mean_of_mean;
			d_input[t*nchans + c] = MEAN;
		}
	}

	int c = blockIdx.x * blockDim.x + threadIdx.x;

	if(c < nchans) {

		count = 0;
		iters = 0;

		stdev = 1000000.0f;
		mean = MEAN;
		mean_last = 0.0f;
		cutoff = (CUT * stdev); 

		while(abs(mean - mean_last) > ACC) {
			sum = 0.0f;
			sum_squares = 0.0f;
			count = 0;

			for(int t = 0; t < nsamp; t++) {
				float data=(float)d_input[t*nchans + c];
				if(data < (mean + cutoff) && data > (mean - cutoff) ) {
					sum+=data;
					sum_squares+=(data*data);
					count++;
				}
			}
			mean_last = mean;
			mean = (sum/(float)count);
			sum_squares = ((sum_squares / count) - (mean * mean));
		    	stdev = sqrt(sum_squares);
			cutoff = (CUT * stdev); 

			iters++;
			if(iters > ITER) break;
		}

		for(int t = 0; t < nsamp; t++) {
			if((d_input[t*nchans + c]-mean < R_CUT*stdev) && (d_input[t*nchans + c]-mean > - R_CUT*stdev)) {
				d_input[t*nchans + c] = (unsigned short)((float)d_input[t*nchans + c]-(float)mean+MEAN);
			} else {
				d_input[t*nchans + c] = MEAN;
			}
		}
	}
}
//}}}

#endif

