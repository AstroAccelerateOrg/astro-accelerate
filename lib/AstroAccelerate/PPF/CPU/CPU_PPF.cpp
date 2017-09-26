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

#include <omp.h>
#include <xmmintrin.h>
#include <immintrin.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mkl.h>

#include <sys/time.h>
#include <unistd.h>

#include "params.h"



#define FpR 8
#define RpCl 2


bool check_errors=0;

double elapsedTime (void)
{
    struct timeval t;
    gettimeofday(&t, 0);
    return ((double)t.tv_sec + ((double)t.tv_usec / 1000000.0));
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


double FIR_check_uni(float *input_data, float *spectra_GPU, float *coeff, int nTaps, int nChannels, int nSpectra){
	nChannels=nChannels*2;
	float ftemp=0;
	float etemp=0;
	double error=0;
	unsigned int count=0;
	for(int bl=0;bl<nSpectra;bl++){
		for(int c=0;c<nChannels;c++){
			ftemp=0;
			for(int t=0;t<nTaps;t++){
				ftemp+=coeff[t*nChannels + c]*input_data[bl*nChannels + t*nChannels + c];
			}//nTaps
			etemp=(ftemp-spectra_GPU[bl*nChannels + c])*(ftemp-spectra_GPU[bl*nChannels + c]);
			error+=etemp;
		}//nChannels
	} //nBlocks

	return(error);
}

// Calculates PPF
void FIR_CPU(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nSpectra){
	int i,c,t,bl,th_id,block_step,remainder,start,end,nThreads,it1,it2;
	nChannels=nChannels*2;
	nSpectra=nSpectra/SpT;
	
	//RpCL = Registers per Cacheline
	//FpR  = Floats per Register
	//CpB = Cachelines per block
	//SpT = Spectra per Thread
	int REG=nChannels/(RpCl*FpR);
	__m256 i_data[RpCl];// because of the size of the cacheline 64byte
	__m256 i_coeff[RpCl];
	__m256 i_spectra[SpT*RpCl];
	__m256 i_temp[RpCl];
	
	int OUTER=(int) REG/CpB;//nuof repetitions
	int REMAINDER,BLOCKS;
	if (OUTER*CpB==REG) {REMAINDER=0;BLOCKS=OUTER;}
	else {REMAINDER=REG-OUTER*CpB;BLOCKS=(OUTER+1);}

	#pragma omp parallel shared(input_data,spectra,coeff) private(it1,it2,OUTER,i,start,end,th_id,bl,t,c,i_data,i_coeff,i_spectra,i_temp)
	{
		th_id = omp_get_thread_num();	// Thread Id
		nThreads = omp_get_num_threads();	// Number of threads
		block_step=nSpectra/nThreads;	// size of the block which thread operates on
		remainder=nSpectra-block_step*nThreads;	// in case nSpectra in not a multiple of nThreads
		
		start=th_id*block_step+remainder;		// Here we calculate start and end of the thread's block
		end=(th_id+1)*block_step+remainder;
		if(th_id==(nThreads-1)) end=nSpectra;
		
		if(th_id<remainder){	// if there is some remainder it is redistributed among threads.
			start=th_id*(block_step+1);
			end=(th_id+1)*(block_step+1);
			if(th_id==(nThreads-1)) end=nSpectra;
		}
		
		for(i=0;i<BLOCKS;i++){ // outer loop through spectra (S_B in article)
			OUTER=CpB;
			if(i==(BLOCKS-1) && REMAINDER!=0) OUTER=REMAINDER;
			for(bl=start;bl<end;bl++){
				for(c=0;c<OUTER;c++){ // outer loop through channels (C_B in article)
					#pragma unroll
					for(it1=0;it1<SpT*RpCl;it1++){
						i_spectra[it1]=_mm256_setzero_ps();
					}
					
					for(t=0;t<nTaps;t++){ // order of the loops is reversed to inner loops (through S_B and C_B) are within tap loop
						#pragma unroll
						for(it2=0;it2<RpCl;it2++){
							i_coeff[it2]=_mm256_load_ps(&coeff[(RpCl*c+it2)*FpR + t*nChannels + i*CpB*RpCl*FpR]);
						}
						
						#pragma unroll
						for(it1=0;it1<SpT;it1++){ // Inner loops for C_B an d S_B
							#pragma unroll
							for(it2=0;it2<RpCl;it2++){
								i_data[it2]=_mm256_load_ps(&input_data[(t+SpT*bl+it1)*nChannels + (RpCl*c+it2)*FpR + i*CpB*RpCl*FpR]);
								i_temp[it2]=_mm256_mul_ps(i_data[it2],i_coeff[it2]);
								i_spectra[(it1*RpCl)+it2]=_mm256_add_ps(i_temp[it2],i_spectra[(it1*RpCl)+it2]);
							}
						}
					} // for nTaps

					#pragma unroll
					for(it2=0;it2<RpCl;it2++){
						#pragma unroll
						for(it1=0;it1<SpT;it1++){
							_mm256_stream_ps(&spectra[(RpCl*c+it2)*FpR + (SpT*bl+it1)*nChannels + i*CpB*RpCl*FpR],i_spectra[(it1*RpCl)+it2]);
						}
					}
					
				}// nChannels
			}// for nBlocks
		} // for BLOCKS
	} // parallel block
}

// This returns Hanning window coefficient.
float Hanning_window(int t, int nChannels){
	return(0.54-0.46*cos(2.0*3.141592654*t/(nChannels-1)));
}

// PPF coefficients used in the article. Coefficients are a sinc(x) function with Hanning window applied to it.
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

int main(int argc, char* argv[])
{
	kmp_set_defaults("KMP_AFFINITY=compact");
	omp_set_num_threads(NUMTHREADS);
	//******************************* ARGUMENTS ********************************
	if (argc!=4) {
		printf("Argument error!\n");
		printf(" nTaps nChannels nSpectra \n");
        return 1;
    }
	char * pEnd;
	
	int nTaps=strtol(argv[1],&pEnd,10);
	int nChannels=strtol(argv[2],&pEnd,10);
	int nSpectra=strtol(argv[3],&pEnd,10);
	const int nRuns=10;
	int nSpectra_to_run;
	int itemp;
	
	if (nSpectra%SpT!=0) {
		itemp=(int) (nSpectra/SpT);
		nSpectra_to_run=(itemp+1)*SpT;
	}
	else nSpectra_to_run=nSpectra;
	//******************************* ARGUMENTS ********************************
	
	int input_size=(nSpectra_to_run+nTaps-1)*nChannels*2;
	int output_size=nSpectra_to_run*nChannels*2;
	int coeff_size=nTaps*nChannels*2;

	float *input;
	float *coeff;
	float *output;
	
	//*************************** DATA GENERATION ******************************
	int data_size,f;
	data_size=(nSpectra_to_run+nTaps-1)*nChannels*2;
	
	input = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	output = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	coeff = (float*)_mm_malloc( nChannels * nTaps * 2 * sizeof(float) ,64);
	
	srand(time(NULL));
	for(f=0;f<data_size;f++){
		input[f]=rand() / (float)RAND_MAX;
	}
	
	PPF_coefficients(coeff,nChannels,nTaps,nChannels);
	//*************************** DATA GENERATION ******************************
	
	
	// -------> Time measuring
	double FIR_start, FIR_time, FIR_sum, FIR_mean;
	double FFT_start, FFT_time, FFT_sum, FFT_mean;
	double Total_time;
    // -------> Time measuring
	
	int count,r;
	bool debug=true;
	
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	FIR_sum=0; FFT_sum=0;
	for(r=0;r<nRuns;r++){
		// FIR
		FIR_start=omp_get_wtime();
		FIR_CPU(input,output,coeff,nTaps,nChannels,nSpectra_to_run);
		FIR_time = omp_get_wtime() - FIR_start;
		
		if (check_errors) printf("Error: %e\n",FIR_check_uni(input,output,coeff,nTaps,nChannels,nSpectra));
		
		//FFT
		FFT_start=omp_get_wtime();
		FFT_separated_single((MKL_Complex8 *) output,nSpectra,nChannels);
		FFT_time = omp_get_wtime() - FFT_start;
		
		FIR_sum+=FIR_time;
		FFT_sum+=FFT_time;
	}
	FIR_mean=FIR_sum/nRuns;
	FFT_mean=FFT_sum/nRuns;
	
	printf("%d %d %d %0.6f %0.6f\n",nSpectra,nChannels,nTaps,FIR_mean,FFT_mean);
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	
	
	
	//******************************* CLEAN UP *********************************
	_mm_free(input);
	_mm_free(coeff);
	_mm_free(output);
	//******************************* CLEAN UP *********************************
    return EXIT_SUCCESS;
}


//Compile:
// icc FIR_filter_autotune_v2.cpp  Filter_window.cpp -o FIR_filter_autotune_v2.exe -openmp -mkl -lfftw3f -include Filter_window.h params.h  -I/home/kadamek/include/ -L/home/kadamek/lib/ -O3
