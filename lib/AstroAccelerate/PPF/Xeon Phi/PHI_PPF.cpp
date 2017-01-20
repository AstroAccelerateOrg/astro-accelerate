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

#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <immintrin.h>
#include <time.h>

#include "params.h"

#define ALIGN_SIZE   (64)


int check_errors=0;

double FIR_check_uni(float *input_data, float *spectra_GPU, float *coeff, int nTaps, int nChannels, int nSpectra){
	nChannels=nChannels*2;
	double ftemp=0;
	double etemp=0;
	double error=0;
	for(int bl=0;bl<nSpectra;bl++){
		for(int c=0;c<nChannels;c++){
			ftemp=0;
			for(int t=0;t<nTaps;t++){
				ftemp+=coeff[t*nChannels + c]*input_data[bl*nChannels + t*nChannels + c];
			}//nTaps
			etemp=abs(ftemp-spectra_GPU[bl*nChannels + c]);
			error+=etemp;
		}//nChannels
	} //nBlocks

	return(error);
}


void FFT_separated_single(MKL_Complex8 *input, int nSpectra, int nChannels){
	MKL_LONG status;
	int nThreads;
	int th=0;

	DFTI_DESCRIPTOR_HANDLE handle = 0;
	status=DftiCreateDescriptor(&handle, DFTI_SINGLE, DFTI_COMPLEX,1,(MKL_LONG) nChannels);
	nThreads = omp_get_max_threads ();
	status=DftiSetValue(handle, DFTI_PLACEMENT, DFTI_INPLACE);
	status=DftiSetValue (handle, DFTI_NUMBER_OF_USER_THREADS, nThreads);
	status=DftiCommitDescriptor(handle);
	
	#pragma omp parallel for private(th,status) num_threads(nThreads)
	for (th = 0; th < nSpectra; th++) {
		status=DftiComputeForward(handle, &input[th*nChannels]);
		if (status && !DftiErrorClass(status,DFTI_NO_ERROR)) {
		   printf ("Error: %s\n", DftiErrorMessage(status));
		}
	}
		
	DftiFreeDescriptor(&handle);
}

//nSpectra=SpT*TpC*nCores 120640 nebo 60320 nebo 30160
void FIR_PHI(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nSpectra){
	int c,i,o,t,bl,th_id,block_step,remainder,nThreads,blt,nCores,Core_start,Core_end,co_id,it1,it2;
	nChannels=nChannels*2;
	nSpectra=nSpectra/SpT;//pre-division for multiple spectra per thread
	nSpectra=nSpectra/4; //pre-division per core
	

	__m512 i_data[SpT*RpCl];
	__m512 i_coeff[RpCl];
	__m512 i_spectra[SpT*RpCl];

	//Assumptions:
	// nSpectra are multiple of SpT*TpC*nCores
	// nChannels are multiple of CpB*RpCl
	
	int REMAINDER,OUTER,INNER,REG;
	REG=nChannels/(RpCl*FpR);
	REMAINDER=nChannels-REG*RpCl*FpR;
	if(REMAINDER!=0){printf("Number of channels is not divisible by 8!\n");exit(101);}
	if(CpB>REG){printf("CpB is set too high\n");exit(103);}
	OUTER=REG/CpB;
	REMAINDER=REG-OUTER*CpB;
	if(REMAINDER!=0) OUTER++;

	#pragma omp parallel shared(input_data,spectra,coeff) private(co_id,Core_start,Core_end,INNER,i,o,blt,th_id,block_step,bl,t,c,i_data,i_coeff,i_spectra,it1,it2)
	{
		th_id = omp_get_thread_num();	// Thread Id
		co_id = (int) th_id/4;			// Core Id assumes compact thread alignment
		nThreads = omp_get_num_threads();
		nCores=(int) nThreads/4;
		block_step=nSpectra/nCores;
		remainder=nSpectra-block_step*nCores;
		Core_start=co_id*block_step;	// Blocks of data are assigned to cores not to threads as in CPU code
		Core_end=(co_id+1)*block_step;
		if(co_id==(nCores-1)) Core_end=nSpectra; 

		th_id=th_id-co_id*4;
		for(o=0;o<OUTER;o++){	// outer loop through channels
			INNER=CpB;
			if(o==OUTER-1 && REMAINDER>0) INNER=REMAINDER;
			for(blt=0;blt<block_step;blt++){ // outer loop through spectra
				bl=Core_start*4 + blt*4 + th_id;
				for(i=0;i<INNER;i++){
					c=o*CpB + i;
					
					#pragma ivdep
					#pragma unroll
					for(it1=0;it1<SpT*RpCl;it1++){
						i_spectra[it1]=_mm512_setzero_ps();
					}
					
					for(t=0;t<nTaps;t++){	// loop through taps. Loops through C_B and S_B are inside this loop
						
						#pragma ivdep
						#pragma unroll
						for(it2=0;it2<RpCl;it2++){
							i_coeff[it2]=_mm512_load_ps(&coeff[(RpCl*c+it2)*FpR+t*nChannels]);
						}
						
						#pragma simd
						#pragma ivdep
						#pragma unroll
						for(it1=0;it1<SpT;it1++){
							#pragma simd
							#pragma ivdep
							#pragma unroll
							for(it2=0;it2<RpCl;it2++){
								i_data[(it1*RpCl)+it2]=_mm512_load_ps(&input_data[(t+SpT*bl+it1)*nChannels+(RpCl*c+it2)*FpR]);
								i_spectra[(it1*RpCl)+it2]=_mm512_fmadd_ps(i_coeff[it2],i_data[(it1*RpCl)+it2],i_spectra[(it1*RpCl)+it2]);
							}
						}
						
					} // for nTaps
					#pragma simd
					#pragma ivdep
					#pragma unroll
					for(it1=0;it1<SpT;it1++){
						#pragma simd
						#pragma ivdep
						#pragma unroll
						for(it2=0;it2<RpCl;it2++){
							_mm512_store_ps(&spectra[(SpT*bl+it1)*nChannels+(RpCl*c+it2)*FpR],i_spectra[(it1*RpCl)+it2]);
						}
					}
					//_mm512_store_ps(&spectra[(RpCl*c)*FpR + bl*nChannels],i_spectra[0]);
				}// for CpB
			}// for nSpectra
		}// for OUTER
	}
	
}

int main(int argc, char* argv[]) {
	//******************************* ARGUMENTS ********************************
	if (argc!=4) {
		printf("Argument error!\n");
		printf("Parameters: nSpectra nChannels nTaps\n");
        return 1;
    }
	char * pEnd;
	
	int nTaps=strtol(argv[3],&pEnd,10);
	int nChannels=strtol(argv[2],&pEnd,10);
	int nSpectra=strtol(argv[1],&pEnd,10);
	const int nRuns=10;
	int nSpectra_to_run;
	
	if (nSpectra%SpT!=0) {
		itemp=(int) (nSpectra/SpT);
		nSpectra_to_run=(itemp+1)*SpT;
	}
	else nSpectra_to_run=nSpectra;
	//******************************* ARGUMENTS ********************************	
	
	//********************************* ALLOCATIONS **********************************
	int bl,t,r;
	int input_size=(nSpectra_to_run+nTaps-1)*nChannels*2;
	int output_size=nSpectra_to_run*nChannels*2;
	int coeff_size=nTaps*nChannels*2;
	
	float  *input;
	float  *output;
	float  *coeff;
	
	input = (float*)_mm_malloc( input_size*sizeof(float) ,64);
	output = (float*)_mm_malloc( output_size*sizeof(float) , 64);
	coeff = (float*)_mm_malloc( coeff_size * sizeof(float) , 64);
	//********************************* ALLOCATIONS **********************************

	//*************************** DATA GENERATION ******************************
	srand (time(NULL));
	// initialize
	#pragma omp parallel for
	for(t=0;t<input_size;t++) {
		input[f]=rand() / (float)RAND_MAX;
	}
	#pragma omp parallel for
	for(t=0;t<output_size;t++) {
		output[t]=0;
	}
	
	PPF_coefficients(coeff,nChannels,nTaps,nChannels);
	//*************************** DATA GENERATION ******************************
	
	// -------> Time measuring
	double FIR_start, FIR_time;
	double FFT_start, FFT_time;
	double Total_time;
    // -------> Time measuring
	
	//********************************************************
	//*                     Performance
	//********************************************************
	
	//********************************************************
	//************  WARM UP
	//*****
	FIR_MK9(input, output, coeff, nTaps, nChannels, nSpectra_to_run);
	FFT_separated_single((MKL_Complex8 *) output,nSpectra,nChannels);
	
	
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	FIR_sum=0;
	FFT_sum=0;
	for(r=0;r<nRuns;r++){
		FIR_start=omp_get_wtime();
		FIR_PHI(input, output, coeff, nTaps, nChannels, nSpectra_to_run);
		FIR_time = omp_get_wtime() - FIR_start;
		
		if (check_errors) printf("Error: %e\n",FIR_check_uni(input,output,coeff,nTaps,nChannels,nSpectra));
		
		FFT_start=omp_get_wtime();
		FFT_separated_single((MKL_Complex8 *) output,nSpectra,nChannels);
		FFT_time = omp_get_wtime() - FFT_start;
		
		FIR_sum+=FIR_time;
		FFT_sum+=FFT_time;
	}
	FIR_mean=FIR_sum/nRuns;
	FFT_mean=FFT_sum/nRuns;
	
	printf("%d %d %d %0.6f %0.6f\n",nSpectra,nTaps,nChannels,FIR_mean,FFT_mean);
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	
	//******************************* CLEAN UP *********************************
	_mm_free(input);
	_mm_free(output);
	_mm_free(coeff);
	//******************************* CLEAN UP *********************************

  return 0;
}
