#ifndef ZERODM_OUTLIERS_KERNEL_H_
#define ZERODM_OUTLIERS_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"


//{{{ zero dm kernel - needs cleaning and optimizing // WA 21/10/16
__global__ void zero_dm_outliers_kernel(unsigned short *d_input, int nchans, int nsamp)
{

	int t  = blockIdx.x * blockDim.x + threadIdx.x;

	int count =0;

	float stdev = 1000000.0f;
	float mean = 0.0f;
	float sum = 0.0f;
	float sum_squares = 0.0f;
	float cutoff = (3.0f * stdev); 

	for(int out=0; out<4; out++) {
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
		mean = (sum/(float)count);
		sum_squares = ((sum_squares / count) - (mean * mean));
	    	stdev = sqrt(sum_squares);
		cutoff = (3.0f * stdev); 
	}

/*	if(mean > 5.0f * stdev) {
//	if(mean > 125.0f && mean < 130.0f) {
		for(int c = 0; c < nchans; c++) {
			float data=(float)d_input[t*nchans + c];

			if(data < (mean - cutoff) || data > (mean + cutoff)) {
				d_input[t*nchans + c]=(unsigned short)128;
			} else {
				d_input[t*nchans + c]=(unsigned short)((unsigned char)((float)d_input[t*nchans + c]-mean+128));
			}
		}
	} else {
		for(int c = 0; c < nchans; c++) {
			d_input[t*nchans + c]=(unsigned short)(128);
		}
	}
*/
//	if(mean > 5.0f * stdev) {
	if(mean > 125.0f && mean < 130.0f) {
		for(int c = 0; c < nchans-4; c++) {
			//float data=(float)d_input[t*nchans + c];
			float data=0.0f;
			for(int x = 0; x<4; x++) data+=(float)d_input[t*nchans + c + x];
			data=data*0.25f;

			if(data < (mean - cutoff) || data > (mean + cutoff)) {
				d_input[t*nchans + c]=(unsigned short)0;
			} else {
				d_input[t*nchans + c]=(unsigned short)((unsigned char)((float)d_input[t*nchans + c]-mean));
			}
		}
	} else {
		for(int c = 0; c < nchans; c++) {
			d_input[t*nchans + c]=(unsigned short)(0);
		}
	}

/*
	if(mean > 5.0f * stdev) {
//	if(mean > 125.0f && mean < 130.0f) {
		for(int c = 0; c < nchans; c++) {
			float data=(float)d_input[t*nchans + c];

			if(data < (mean - cutoff) || data > (mean + cutoff)) {
				d_input[t*nchans + c]=(unsigned short)0;
			} else {
				d_input[t*nchans + c]=(unsigned short)((unsigned char)(((float)d_input[t*nchans + c]-mean)/stdev));
			}
		}
	} else {
		for(int c = 0; c < nchans; c++) {
			d_input[t*nchans + c]=(unsigned short)(0);
		}
	}
*/
/*	if(mean > 5.0f * stdev) {
//	if(mean > 125.0f && mean < 130.0f) {
		for(int c = 0; c < nchans; c++) {
			float data=(float)d_input[t*nchans + c];

			if(data < (mean - cutoff) || data > (mean + cutoff)) {
				d_input[t*nchans + c]=(unsigned short)(128.0f/stdev);
			} else {
				d_input[t*nchans + c]=(unsigned short)((unsigned char)(((float)d_input[t*nchans + c]-mean+128)/stdev));
			}
		}
	} else {
		for(int c = 0; c < nchans; c++) {
			d_input[t*nchans + c]=(unsigned short)(128.0f/stdev);
		}
	}
*/
}
//}}}

#endif

