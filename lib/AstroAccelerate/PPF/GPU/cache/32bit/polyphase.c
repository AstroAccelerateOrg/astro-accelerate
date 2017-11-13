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
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

void Normalize_FFT(float2 *output, int nChannels, int nSpectra, double factor);

void Perform_FFT(float2 *spectra, int nChannels);

void FIR_check(float2 *input_data, float2 *spectra_GPU, float *coeff, int nTaps, int nChannels, int nSpectra, double *cumulative_error, double *mean_error);

void FIR_FFT_check(float2 *input_data, float2 *spectra_GPU, float *coeff, int nTaps, int nChannels, int nSpectra, double *cumulative_error, double *mean_error);

void GPU_Polyphase(float2 *input, float2 *output, float *coeff, int nChannels, int nTaps, int nSpectra);

int Max_columns_in_memory(int nTaps, int nChannels);

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

	int nColumns=Max_columns_in_memory(nTaps,nChannels);
	if (nSpectra < nColumns ) nColumns = nSpectra;
	
	int input_size=(nColumns+nTaps-1)*nChannels;
	int output_size=nColumns*nChannels;
	int coeff_size=nTaps*nChannels;
	
	double cumulative_error,mean_error;

	if (DEBUG) printf("\t\tWelcome\n");

	float2 *h_input;
	float2 *h_output;
	float *h_coeff;

	if (DEBUG) printf("\nHost memory allocation...\t");
		h_input 	= (float2 *)malloc(input_size*sizeof(float2));
		h_output 	= (float2 *)malloc(output_size*sizeof(float2));
		h_coeff 	= (float *)malloc(coeff_size*sizeof(float));
	if (DEBUG) printf("done.");

	if (DEBUG) printf("\nHost memory memset...\t\t");
		memset(h_output, 0.0, output_size*sizeof(float2));	
	if (DEBUG) printf("done.");

	if (DEBUG) printf("\nLoad window coefficients...\t");
		//Load_window_data(h_coeff);
		srand(time(NULL));
		for (int i = 0; i < coeff_size; i++)
			h_coeff[i] = rand() / (float)RAND_MAX;
	if (DEBUG) printf("done.");


	if (DEBUG) printf("\nRandom data set...\t\t");	
		srand(time(NULL));
		for (int i=0; i < (int) input_size; i++){
			h_input[i].x = rand() / (float)RAND_MAX;
			h_input[i].y = rand() / (float)RAND_MAX;
		}
	if (DEBUG) printf("done.");

	GPU_Polyphase(h_input, h_output, h_coeff, nChannels, nTaps, nSpectra);
	/*
	if (DEBUG && CHECK){ 
		printf("\nTesting FIR...\t\t");
		FIR_check(h_input, h_output, h_coeff, nTaps, nChannels, nSpectra, &cumulative_error, &mean_error);
		printf("Cumulative Error: %e, Mean Error: %e;",cumulative_error,mean_error);
	}
	*/
	if (DEBUG && CHECK){
		printf("\nTesting FIR+FFT...\t");
		FIR_FFT_check(h_input, h_output, h_coeff, nTaps, nChannels, nSpectra, &cumulative_error, &mean_error);
		printf("Cumulative Error: %e, Mean Error: %e;",cumulative_error,mean_error);
	}

	delete[] h_input;
	delete[] h_output;
	delete[] h_coeff;
	
	cudaDeviceReset();
	
	if (CHECK && mean_error>1.0e-4) exit(1);

	if (DEBUG) printf("\nThat's All folks!\n");

	return (0);
}
