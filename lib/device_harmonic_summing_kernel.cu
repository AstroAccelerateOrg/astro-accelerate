

#ifndef HARMONIC_SUMMING_KERNEL_H_
#define HARMONIC_SUMMING_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"


__global__ void PHS_GPU_kernel_old(float const* __restrict__ d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics){
	float signal_mean=d_MSD[0];
	float signal_sd=d_MSD[1];
	float HS_value, temp_SNR, SNR;
	ushort max_SNR_harmonic;
	int pos;

	// reading 0th harmonic, i.e. fundamental frequency
	pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
	if( (blockIdx.y*blockDim.x + threadIdx.x)<nSpectra ){
		HS_value = __ldg(&d_input[pos]);
		SNR = (HS_value - signal_mean)/(signal_sd);
		max_SNR_harmonic = 0;
		
		if(blockIdx.x>0) {
			for(int f=1; f<nHarmonics; f++){
				if( (blockIdx.x + f*blockIdx.x)<nTimesamples ){
					pos = (blockIdx.x + f*blockIdx.x)*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
					HS_value = HS_value + __ldg(&d_input[pos]);
					temp_SNR = __frsqrt_rn(f+1)*(HS_value - signal_mean*(f+1))/(signal_sd); //assuming white noise 
					if(temp_SNR > SNR){
						SNR = temp_SNR;
						max_SNR_harmonic = f;
					}
				}
			}
		}
		
		pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
		d_output_SNR[pos] = SNR;
		d_output_harmonics[pos] = max_SNR_harmonic;
	}
}


__global__ void PHS_GPU_kernel(float const* __restrict__ d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics){
	//extern __shared__ float s_MSD[];
	float HS_value, temp_SNR, SNR;
	ushort max_SNR_harmonic;
	int pos;

	//int itemp = (2*nHarmonics)/blockDim.x;
	//for(int f=0; f<itemp; f++){
	//	pos = f*blockDim.x + threadIdx.x;
	//	if(pos<nHarmonics) s_MSD[pos] = d_MSD[pos];
	//}

	// reading 0th harmonic, i.e. fundamental frequency
	pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
	if( (blockIdx.y*blockDim.x + threadIdx.x)<nSpectra ){
		HS_value = __ldg(&d_input[pos]);
		SNR = (HS_value - __ldg(&d_MSD[0]))/(__ldg(&d_MSD[1]));
		max_SNR_harmonic = 0;
		
		if(blockIdx.x>0) {
			for(int f=1; f<nHarmonics; f++) {
				if( (blockIdx.x + f*blockIdx.x)<nTimesamples ) {
					pos = (blockIdx.x + f*blockIdx.x)*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
					HS_value = HS_value + __ldg(&d_input[pos]);
					temp_SNR = (HS_value - __ldg(&d_MSD[f*2]))/(__ldg(&d_MSD[2*f+1]));
					if(temp_SNR > SNR) {
						SNR = temp_SNR;
						max_SNR_harmonic = f;
					}
				}
			}
		}
		
		pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
		d_output_SNR[pos] = SNR;
		d_output_harmonics[pos] = max_SNR_harmonic;
	}
}

#endif