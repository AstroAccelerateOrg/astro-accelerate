/*************************************************************************
    This is GPU implementation of a polyphase filter. 
    Copyright (C) 2015  Adamek Karel, Novotny Jan, Armour Wes

    This file is part of Astro-Accelerate PolyPhase Filter (AAPPF).

    AAPPF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AAPPF is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AAPPF.  If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/

#include "debug.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

void Normalize_FFT(float2 *output, int nChannels, int nSpectra, double factor);

void Perform_FFT(float2 *spectra, int nChannels);

void FIR_check(uchar2 *input_data, float2 *spectra_GPU, float *coeff, int nTaps, int nChannels, int nSpectra, double *cumulative_error, double *mean_error);

void FIR_FFT_check(uchar2 *input_data, float2 *spectra_GPU, float *coeff, int nTaps, int nChannels, int nSpectra, double *cumulative_error, double *mean_error);

void GPU_Polyphase(uchar2 *input, float2 *output, float *coeff, int nChannels, int nTaps, int nSpectra);

int Max_columns_in_memory(int nTaps, int nChannels);

float Hanning_window(int t, int nChannels){
	return(0.54-0.46*cos(2.0*3.141592654*t/(nChannels-1)));
}

void PPF_coefficients(float *coeff, int nChannels, int nTaps, int width){
	float mag,x,ftemp;
	int f,start,end,zero;
	zero=nTaps*nChannels/2;
	start=zero-width;end=zero+width;
	for(f=0;f<nChannels*nTaps;f++) {
		x=(float) f-zero;
		ftemp=sin(x*3.141592654/width)/(x*3.141592654/width);
		if(ftemp==ftemp) mag=ftemp;
		mag=mag*Hanning_window(f,nTaps*nChannels);
		coeff[f]=mag;
	}
}
int main(int argc, char* argv[]) {
	
	if (argc!=4) {
		printf("Argument error!\n");
		printf(" nChannels nTaps nSpectra \n");
        return 1;
    }
	char * pEnd;
	
	int nChannels=strtol(argv[1],&pEnd,10);
	int nTaps=strtol(argv[2],&pEnd,10);
	int nSpectra=strtol(argv[3],&pEnd,10);

	int input_size=(nSpectra+nTaps-1)*nChannels;
	int output_size=nSpectra*nChannels;
	int coeff_size=nTaps*nChannels;
	
	double cumulative_error,mean_error;

	if (DEBUG) printf("\t\tWelcome\n");

	uchar2 *h_input;
	float2 *h_output;
	float *h_coeff;

	if (DEBUG) printf("\nHost memory allocation...\t");
		h_input 	= (uchar2 *)malloc(input_size*sizeof(uchar2));
		h_output 	= (float2 *)malloc(output_size*sizeof(float2));
		h_coeff 	= (float *)malloc(coeff_size*sizeof(float));
	if (DEBUG) printf("done.");

	if (DEBUG) printf("\nHost memory memset...\t\t");
		memset(h_output, 0.0, output_size*sizeof(float2));	
	if (DEBUG) printf("done.");

	if (DEBUG) printf("\nLoad window coefficients...\t");
		PPF_coefficients(h_coeff,nChannels,nTaps,nChannels);
	if (DEBUG) printf("done.");


	if (DEBUG) printf("\nRandom data set...\t\t");	
		srand(time(NULL));
		for (int i=0; i < (int) input_size; i++){
			h_input[i].x = rand() % 255;
			h_input[i].y = rand() % 255;
		}
	if (DEBUG) printf("done.");

	GPU_Polyphase(h_input, h_output, h_coeff, nChannels, nTaps, nSpectra);
	/*
	if (DEBUG && CHECK){ 
		printf("\nTesting FIR...\t\t");
		FIR_check((uchar2*) h_input, h_output, h_coeff, nTaps, nChannels, nSpectra, &cumulative_error, &mean_error);
		printf("Cumulative Error: %e, Mean Error: %e;",cumulative_error,mean_error);
	}
	*/
	if (DEBUG && CHECK){
		printf("\nTesting FIR+FFT...\t");
		FIR_FFT_check((uchar2*) h_input, h_output, h_coeff, nTaps, nChannels, nSpectra, &cumulative_error, &mean_error);
		printf("Cumulative Error: %e, Mean Error: %e;",cumulative_error,mean_error);
	}

	delete[] h_input;
	delete[] h_output;
	delete[] h_coeff;
	
	cudaDeviceReset();
	

	if (DEBUG) printf("\nThat's All folks!\n");

	return (0);
}
