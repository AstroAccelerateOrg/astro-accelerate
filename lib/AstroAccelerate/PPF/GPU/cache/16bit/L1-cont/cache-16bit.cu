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

#include <cufft.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include "../timer.h"
#include "../utils_cuda.h"
#include "../utils_file.h"

#include "params.h"
#include "../debug.h"

int device=0;

// --------> maxwell and kepler
__global__ void Fir_L1(short2 const* __restrict__ d_data, float2* d_spectra, float const* __restrict__ d_coeff, int nTaps, int nChannels, int yshift) {
	int t = 0;
	int bl= SPECTRA_PER_BLOCK*blockIdx.x*nChannels;
	float temp;
	short2 utemp;
	float2 ftemp[SPECTRA_PER_BLOCK];
	
	for (int i = 0; i<SPECTRA_PER_BLOCK;i++){
	  ftemp[i].x = 0.0f;
	  ftemp[i].y = 0.0f;
	}

	for(t=threadIdx.x + yshift; t<(nTaps)*nChannels; t+=nChannels){
	  temp = __ldg(&d_coeff[t]); 
	  for(int i=0;i<SPECTRA_PER_BLOCK;i++){
	    utemp = __ldg(&d_data[bl+i*nChannels + t]);
	    ftemp[i].x += temp*((float)(utemp.x));
	    ftemp[i].y += temp*((float)(utemp.y));
	  }
	}

	t=bl + threadIdx.x + yshift;
	for (int i=0;i<SPECTRA_PER_BLOCK;i++){
		d_spectra[t + i*nChannels]=ftemp[i];
	}

	return;
}

// --------------> fermi generation
/*
__global__ void Fir_L1(short2 const* __restrict__ d_data, float2* d_spectra, float const* __restrict__ d_coeff, int nTaps, int nChannels, int yshift) {
	int t = 0;
	int bl= SPECTRA_PER_BLOCK*blockIdx.x*nChannels;
	float temp;
	float2 ftemp[SPECTRA_PER_BLOCK];
	
	for (int i = 0; i<SPECTRA_PER_BLOCK;i++){
	  ftemp[i].x = 0.0f;
	  ftemp[i].y = 0.0f;
	}

	for(t=threadIdx.x + yshift; t<(nTaps)*nChannels; t+=nChannels){
  	  temp = d_coeff[t]; 
	  for(int i=0;i<SPECTRA_PER_BLOCK;i++){
	    ftemp[i].x += temp*((float)d_data[bl+i*nChannels + t].x);
	    ftemp[i].y += temp*((float)d_data[bl+i*nChannels + t].y);
	  }
	}

	t=bl + threadIdx.x + yshift;
	for (int i=0;i<SPECTRA_PER_BLOCK;i++){
		d_spectra[t + i*nChannels]=ftemp[i];
	}

	return;
}
*/

int Check_parameters(int nTaps, int nChannels, int nSpectra) {
	long int nColumns, maxgrid_x,itemp;
	size_t free_mem,total_mem;
	cudaDeviceProp devProp;
	
	checkCudaErrors(cudaSetDevice(device));
	checkCudaErrors(cudaGetDeviceProperties(&devProp,device));
	maxgrid_x = devProp.maxGridSize[0];
	cudaMemGetInfo(&free_mem,&total_mem);
	
	nColumns=(free_mem-nChannels*sizeof(float)*nTaps -(nTaps-1)*nChannels*sizeof(float2))/(2.0*sizeof(float2)*nChannels);
	if(maxgrid_x*SPECTRA_PER_BLOCK<nColumns) nColumns=maxgrid_x*SPECTRA_PER_BLOCK;
	nColumns=(int) nColumns*0.9;
	itemp=(int) nColumns/SPECTRA_PER_BLOCK;
	nColumns=itemp*SPECTRA_PER_BLOCK;
	return(nColumns);
}


void GPU_Polyphase(short2 *input, float2 *output, float *coeff, int nChannels, int nTaps, int nSpectra){
	int maxgrid_x;
	
	//---------> Initial nVidia stuff
	int devCount;
	cudaDeviceProp devProp;
	size_t free_mem,total_mem;

	checkCudaErrors(cudaGetDeviceCount(&devCount));

	if (DEBUG) {
		printf("\nThere are %d devices.", devCount);
		for (int i = 0; i < devCount; i++){
			checkCudaErrors(cudaGetDeviceProperties(&devProp,i));
			printf("\n\t Using device:\t\t\t%s\n", devProp.name);
			printf("\n\t Max grid size:\t\t\t%d\n", devProp.maxGridSize[1]);
		}
	}

	checkCudaErrors(cudaSetDevice(device));
	checkCudaErrors(cudaGetDeviceProperties(&devProp,device));

	maxgrid_x = devProp.maxGridSize[1];

	cudaMemGetInfo(&free_mem,&total_mem);
	printf("\nDevice has %f MiB of total memory, which %f MiB is available.\n", (double) total_mem/(1024.0*1024.0), (double) free_mem/(1024.0*1024.0));
	
	//---------> Measurements
	double coeff_transport_in=0.0, transfer_in=0.0, transfer_out=0.0, fir_time=0.0, fft_time=0.0;
	GpuTimer timer;
	
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
	
	//---------> Spectra
	long int nColumns;
	int nCUDAblocks,Sremainder,nRepeats,itemp;
	nColumns=(free_mem-nChannels*sizeof(float)*nTaps -(nTaps-1)*nChannels*sizeof(float2))/(2.0*sizeof(float2)*nChannels);
	if(maxgrid_x*SPECTRA_PER_BLOCK<nColumns) nColumns=maxgrid_x*SPECTRA_PER_BLOCK;
	nColumns=(int) nColumns*0.9;
	itemp=(int) nColumns/SPECTRA_PER_BLOCK;
	nColumns=itemp*SPECTRA_PER_BLOCK;
	nRepeats=(int) (nSpectra/nColumns);
	Sremainder=nSpectra-nRepeats*nColumns;
	if (DEBUG) printf("Maximum number of spectra in memory is %d which is %e MB \n",nColumns, (double) (nColumns*nChannels*sizeof(float2)/(1024.0*1024.0))   );
	//if (DEBUG) printf("nSpectra is split into %d chunks. Sremainder: %d\n",nRepeats,Sremainder);
	
	//---------> Channels
	int nKernels=(int) nChannels/THREADS_PER_BLOCK; //Head size
	int Kremainder=nChannels-nKernels*THREADS_PER_BLOCK; //Tail size
	
	//---------> Memory allocation
	if (DEBUG) printf("\nDevice memory allocation...: \t\t");
	int input_size=nChannels*(nColumns+nTaps-1);
	int output_size=nChannels*nColumns;
	int coeff_size=nChannels*nTaps;
	float2 *d_output;
	float *d_coeff;
	short2 *d_input;
	timer.Start();
	checkCudaErrors(cudaMalloc((void **) &d_input,  sizeof(short2)*input_size));
	checkCudaErrors(cudaMalloc((void **) &d_output, sizeof(float2)*output_size));
	checkCudaErrors(cudaMalloc((void **) &d_coeff,  sizeof(float)*coeff_size));
	timer.Stop();
	if (DEBUG) printf("done in %g ms.", timer.Elapsed());
	
	
	//---------> Transfer to the device
	if (DEBUG) printf("\nCopy coefficients from host to device...\t");
	timer.Start();
	checkCudaErrors(cudaMemcpy(d_coeff, coeff, coeff_size*sizeof(float), cudaMemcpyHostToDevice));
	timer.Stop();
	coeff_transport_in=timer.Elapsed();
	if (DEBUG) printf("done in %g ms.\n", coeff_transport_in);
		
	//---------> CUDA block and CUDA grid parameters
	nCUDAblocks=(int) nColumns/SPECTRA_PER_BLOCK;
	dim3 GridSize(nCUDAblocks, 1, 1);
	dim3 BlockSize(THREADS_PER_BLOCK, 1, 1);
	
	//---------> Transfer to the device
	if (DEBUG) printf("\nCopy data from host to device...\t");
	timer.Start();
	if(nRepeats>0){
		checkCudaErrors(cudaMemcpy(d_input, input, input_size*sizeof(short2), cudaMemcpyHostToDevice));
	}
	else {
		checkCudaErrors(cudaMemcpy(d_input, &input[output_size*nRepeats], (Sremainder+nTaps-1)*nChannels*sizeof(short2), cudaMemcpyHostToDevice));
	}
	timer.Stop();
	transfer_in+=timer.Elapsed();
	if (DEBUG) printf("done in %g ms.\n", timer.Elapsed());
		
	// ---------------> Polyphase filter
	for (int r = 0; r < nRepeats; r++){
		//---------> FIR part
		BlockSize.x=THREADS_PER_BLOCK;
		timer.Start();
		for (int nutak=0;nutak<nKernels;nutak++){	
			Fir_L1<<<GridSize, BlockSize>>>((short2*) d_input, (float2*) d_output, d_coeff, nTaps, nChannels, nutak*THREADS_PER_BLOCK);
		}
		if (Kremainder>0){
			BlockSize.x=Kremainder;
			Fir_L1<<<GridSize, BlockSize>>>((short2*) d_input, (float2*) d_output, d_coeff, nTaps, nChannels, nKernels*THREADS_PER_BLOCK);
		}
		timer.Stop();
	 	fir_time+=timer.Elapsed();
	}
	if (Sremainder>0){
		itemp=(int) Sremainder/SPECTRA_PER_BLOCK;
		itemp=(itemp+1)*(SPECTRA_PER_BLOCK);
		nCUDAblocks=itemp/SPECTRA_PER_BLOCK;
		GridSize.x=nCUDAblocks;BlockSize.x=THREADS_PER_BLOCK;
		//---------> FIR
		BlockSize.x=THREADS_PER_BLOCK;
		timer.Start();
		for (int nutak=0;nutak<nKernels;nutak++){	
			Fir_L1<<<GridSize, BlockSize>>>((short2*) d_input, (float2*) d_output, d_coeff, nTaps, nChannels, nutak*THREADS_PER_BLOCK);
		}
		if (Kremainder>0){
			BlockSize.x=Kremainder;
			Fir_L1<<<GridSize, BlockSize>>>((short2*) d_input, (float2*) d_output, d_coeff, nTaps, nChannels, nKernels*THREADS_PER_BLOCK);
		}
		timer.Stop();
		fir_time+=timer.Elapsed();
	}
	
	//----- error check -----
	checkCudaErrors(cudaGetLastError());
	checkCudaErrors(cudaFree(d_input));
	checkCudaErrors(cudaFree(d_coeff));

	//---------> FFT
	for (int r = 0; r < nRepeats; r++){
		cufftHandle plan;
		cufftResult error;
		error = cufftPlan1d(&plan, nChannels, CUFFT_C2C, nColumns);
		if (CUFFT_SUCCESS != error){
			printf("CUFFT error: %d", error);
		}
		
		//execute plan and copy back to host
			timer.Start();
				cufftExecC2C(plan, (cufftComplex *)d_output, (cufftComplex *)d_output, CUFFT_FORWARD);
			timer.Stop();
			fft_time += timer.Elapsed();
		
		//Destroy the cuFFT plan
		cufftDestroy(plan);
	}
	if (Sremainder>0){
		//---------> FFT
		cufftHandle plan;
		cufftResult error;
		error = cufftPlan1d(&plan, nChannels, CUFFT_C2C, Sremainder);
		if (CUFFT_SUCCESS != error){
			printf("CUFFT error: %d", error);
		}
		
		//execute plan and copy back to host
			timer.Start();
				cufftExecC2C(plan, (cufftComplex *)d_output, (cufftComplex *)d_output, CUFFT_FORWARD);
			timer.Stop();
			fft_time += timer.Elapsed();
		
		//Destroy the cuFFT plan
		cufftDestroy(plan);
	}
	
	//---------> Transfer to the host
	if (DEBUG) printf("Copy data from device to host \t");
	timer.Start();
	if(nRepeats>0){
		checkCudaErrors(cudaMemcpy(output,d_output,output_size*sizeof(float2), cudaMemcpyDeviceToHost));
	}
	else {
		checkCudaErrors(cudaMemcpy(&output[output_size*nRepeats],d_output,Sremainder*nChannels*sizeof(float2), cudaMemcpyDeviceToHost));	
	}
	timer.Stop();
	if (DEBUG) printf("done in %g ms.\n", timer.Elapsed());
	transfer_out+=timer.Elapsed();
	
	//----- error check -----
	checkCudaErrors(cudaGetLastError());
	//checkCudaErrors(cudaDeviceSynchronize());
	//-----------------------
	
	if (DEBUG) printf("Number of spectra: %d;\nNumber of Channels: %d;\nNumber of Taps: %d;\nFIR filter execution time: %0.3f ms;\ncuFFT execution time: %0.3f ms;\nPolyphase execution time: %0.3f ms;\nData transfer time %0.3f ms\n",nSpectra,nChannels,nTaps, fir_time, fft_time, fir_time + fft_time, transfer_in + transfer_out);
	
	if (DEBUG && WRITE){ 
		char str[200];
		sprintf(str,"GPU-polyphase.dat");
		printf("\n Write results into file...\t");
		save_time(str, nSpectra, fir_time, fft_time, transfer_in, transfer_out, nChannels, nTaps, SPECTRA_PER_BLOCK, THREADS_PER_BLOCK, 1);
		printf("\t done.\n-------------------------------------\n");
	}
	
	//---------> Cleanup
	if (DEBUG) printf("Free device data.....\t");
	checkCudaErrors(cudaFree(d_output));
	if (DEBUG) printf("\t done.\n-----------------------------------\n");

}
