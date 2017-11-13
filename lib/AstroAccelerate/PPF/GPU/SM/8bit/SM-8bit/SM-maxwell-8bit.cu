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
#include "../debug.h"
#include "../timer.h"
#include "../utils_cuda.h"
#include "../utils_file.h"


#define WARP 32

#include "params.h"

int device=0;
int SM_Columns=(DATA_SIZE/WARP-TAPS+1);

// Maxwell and Kepler
__global__ void Fir_shared_8bit(uchar2 const* __restrict__ d_data, float2* d_spectra, float const* __restrict__ d_coeff, int nChannels) {
	float2 ftemp;
	uchar2 utemp;
	int memblock, localId, s_mempos, g_mempos, num_spectra, start_column, warpId, itemp;
	int tx = threadIdx.x;
	
	__shared__ uchar4 s_data[DATA_SIZE];
	__shared__ float s_coeff[COEFF_SIZE];
	
	warpId=((int) tx/WARP);
	memblock=warpId*SUBBLOCK_SIZE;
	localId= tx - ((int) tx/WARP)*WARP;
	num_spectra = (DATA_SIZE/WARP-TAPS+1);
	

	// read data
	for (int i = 0; i < SUBBLOCK_SIZE; i++){
		start_column = memblock+i;
		if (start_column<DATA_SIZE/WARP) {
			s_mempos = start_column*WARP + localId;
			g_mempos = localId + start_column*nChannels + blockIdx.x*WARP + blockIdx.y*num_spectra*nChannels;
			utemp=__ldg(&d_data[g_mempos]);			
			s_data[s_mempos].x = utemp.x;
			s_data[s_mempos].y = utemp.y;
		}
	}

	itemp=(int) (TAPS/(THREADS_PER_BLOCK/WARP)) + 1;
	for(int f=0;f<itemp;f++){
		start_column=warpId+f*(THREADS_PER_BLOCK/WARP);
		if(start_column<TAPS){
			s_coeff[start_column*WARP+localId]=__ldg(&d_coeff[start_column*nChannels + blockIdx.x*WARP + localId]);
		}
	}
	
	__syncthreads();
	
	for (int i = 0; i < SUBBLOCK_SIZE; i++){
		start_column = memblock+i;
		if (start_column < num_spectra){
			s_mempos = start_column*WARP + localId;
			ftemp.x = 0.0f;ftemp.y = 0.0f;
			for (int j = 0; j < TAPS; j++){
				ftemp.x += s_coeff[j*WARP + localId]*(s_data[s_mempos + j*WARP].x);
				ftemp.y += s_coeff[j*WARP + localId]*(s_data[s_mempos + j*WARP].y);
			} //ntaps
			d_spectra[localId + i*nChannels + memblock*nChannels + blockIdx.x*WARP + blockIdx.y*num_spectra*nChannels] = ftemp;
		} //if
	} // number of columns
	//return;
}

// Fermi
/*
__global__ void Fir_shared_8bit(uchar4 const* __restrict__ d_data, float4* d_spectra, float const* __restrict__ d_coeff, int nChannels) {
	
	float4 ftemp;
	int memblock, localId, s_mempos, g_mempos, num_spectra, start_column, warpId, itemp;
	int tx = threadIdx.x;
	
	__shared__ uchar4 s_data[DATA_SIZE];
	__shared__ float s_coeff[COEFF_SIZE];
		
	warpId=((int) tx/WARP);
	memblock=warpId*SUBBLOCK_SIZE;
	localId= tx - ((int) tx/WARP)*WARP;
	num_spectra = (DATA_SIZE/WARP-TAPS+1);
	

	// read data
	for (int i = 0; i < SUBBLOCK_SIZE; i++){
		start_column = memblock+i;
		if (start_column<DATA_SIZE/WARP) {
			s_mempos = start_column*WARP + localId;
			g_mempos = localId + (blockIdx.y*num_spectra + start_column)*nChannels/2 + blockIdx.x*WARP;
			s_data[s_mempos] = d_data[g_mempos];
		}
	}
	
	itemp=(int) (TAPS/(THREADS_PER_BLOCK/WARP)) + 1;
	for(int f=0;f<itemp;f++){
		start_column=warpId+f*(THREADS_PER_BLOCK/WARP);
		if(start_column<TAPS){
			s_coeff[2*start_column*WARP+localId]=d_coeff[start_column*nChannels + 2*blockIdx.x*WARP + localId];
			s_coeff[2*start_column*WARP+localId+WARP]=d_coeff[start_column*nChannels + (2*blockIdx.x+1)*WARP + localId];
		}
	}
	
	__syncthreads();
	
	for (int i = 0; i < SUBBLOCK_SIZE; i++){
		start_column = memblock+i;
		if (start_column < num_spectra){
			s_mempos = start_column*WARP + localId;
			ftemp.x = 0.0f;ftemp.y = 0.0f;
			ftemp.z = 0.0f;ftemp.w = 0.0f;
			for (int j = 0; j < TAPS; j++){
				ftemp.x += s_coeff[2*j*WARP + 2*localId]*((float) s_data[s_mempos + j*WARP].x);
				ftemp.y += s_coeff[2*j*WARP + 2*localId]*((float) s_data[s_mempos + j*WARP].y);
				ftemp.z += s_coeff[2*j*WARP + 2*localId + 1]*((float) s_data[s_mempos + j*WARP].z);
				ftemp.w += s_coeff[2*j*WARP + 2*localId + 1]*((float) s_data[s_mempos + j*WARP].w);
			} //ntaps
			d_spectra[localId + start_column*nChannels/2 + blockIdx.x*WARP + blockIdx.y*num_spectra*nChannels/2] = ftemp;
		} //if
	} // number of columns
	//return;
}
*/



int Max_columns_in_memory(int nTaps, int nChannels) {
	long int maxColumns, maxgrid_y,itemp;

	size_t free_mem,total_mem;
	cudaDeviceProp devProp;
	
	checkCudaErrors(cudaSetDevice(device));
	checkCudaErrors(cudaGetDeviceProperties(&devProp,device));
	maxgrid_y = devProp.maxGridSize[1];
	cudaMemGetInfo(&free_mem,&total_mem);
	
	maxColumns=((long int) free_mem - nChannels*sizeof(float)*TAPS -(TAPS-1)*nChannels*sizeof(float2))/(2.0*sizeof(float2)*nChannels + sizeof(uchar2)*nChannels);
	maxColumns=(int) (maxColumns*0.9);
	if(maxgrid_y*SM_Columns<maxColumns) maxColumns=maxgrid_y*SM_Columns;
	itemp=(int) (maxColumns/SM_Columns);
	maxColumns=itemp*SM_Columns;
	return(maxColumns);
}

void Polyphase_GPU_init(){
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
}

void Polyphase_GPU_benchmark(uchar2 *d_input, float2 *d_output, float *d_coeff, int nChannels, int nTaps, int nSpectra, double *fir_time, double *fft_time){
	GpuTimer timer;
	//---------> CUDA block and CUDA grid parameters
	int nCUDAblocks_y=(int) (nSpectra/SM_Columns);
	int nCUDAblocks_x=(int) (nChannels/WARP); //Head size
	dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);	//nCUDAblocks_y goes through spectra
	dim3 blockSize(THREADS_PER_BLOCK, 1, 1); 		//nCUDAblocks_x goes through channels
	
	//---------> FIR filter part
	timer.Start();
	Fir_shared_8bit<<<gridSize, blockSize>>>( d_input, d_output, d_coeff, nChannels);
	timer.Stop();
	*fir_time += timer.Elapsed();
	
	//---------> DFT part
	cufftHandle plan;
	cufftResult error;
	error = cufftPlanMany(&plan, 1, &nChannels, &nChannels, 1, nChannels, &nChannels, 1, nChannels, CUFFT_C2C, nSpectra);
	if (CUFFT_SUCCESS != error){
		printf("CUFFT error: %d\n", error);
	}
	
	timer.Start();
	cufftExecC2C(plan, (cufftComplex *)d_output, (cufftComplex *)d_output, CUFFT_FORWARD);
	timer.Stop();
	*fft_time += timer.Elapsed();

	cufftDestroy(plan);
}



void GPU_Polyphase(uchar2 *input, float2 *output, float *coeff, const int nChannels, const int nTaps, const int nSpectra){
	if (nTaps!=TAPS) {printf("nTaps!=TAPS");exit(101);}

	
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
			printf("\n\t Shared mem per block:\t\t%d\n", devProp.sharedMemPerBlock);
		}
	}
	checkCudaErrors(cudaSetDevice(device));
	checkCudaErrors(cudaGetDeviceProperties(&devProp,device));

	
	cudaMemGetInfo(&free_mem,&total_mem);
	if(DEBUG) printf("\nDevice has %ld MB of total memory, which %ld MB is available.\n", (long int) total_mem/(1000*1000), (long int) free_mem/(1000*1000));
	
	//---------> Measurements
	double transfer_in=0.0, transfer_out=0.0, fir_time=0.0, fft_time=0.0;
	GpuTimer timer; // if set before set device getting errors - invalid handle  

	
	//---------> Spectra
	int maxColumns,Sremainder,nRepeats,itemp,Spectra_to_run;
	maxColumns=Max_columns_in_memory(nTaps,nChannels); // Maximum number of columns which fits into memory
	nRepeats=(int) (nSpectra/maxColumns); 
	Sremainder=nSpectra-nRepeats*maxColumns;
	itemp=(int) (Sremainder/SM_Columns);
	if( (Sremainder-itemp*SM_Columns)>0 ) itemp++;
	Spectra_to_run=itemp*SM_Columns; // Since shared memory kernel has fixed number of columns it loads we need to process more spectra then needed.
	//if (DEBUG) printf("Columns per threadblock %d \n",SM_Columns);
	if (DEBUG) printf("Maximum number of spectra in memory is %d which is %e MB \n",maxColumns, (double) (maxColumns*nChannels*sizeof(float2)/(1024.0*1024.0))   );
	//if (DEBUG) printf("nColumns is split into %d chunks of %d spectra and into remainder of %d spectra.\n",nRepeats,maxColumns,Sremainder);
	
	//---------> Channels
	int nCUDAblocks_x=(int) nChannels/WARP; //Head size
	int Cremainder=nChannels-nCUDAblocks_x*WARP; //Tail size
	if(Cremainder>0) {printf("Number of channels must be divisible by 32"); exit(2);}
	
	//---------> Memory allocation
	if (DEBUG) printf("\nDevice memory allocation...: \t\t");
	int input_size=nChannels*(maxColumns+TAPS-1);
	int output_size=nChannels*maxColumns;
	int coeff_size=nChannels*TAPS;

	uchar2 *d_input;
	float2 *d_output;
	float *d_coeff;
	timer.Start();
	checkCudaErrors(cudaMalloc((void **) &d_input,  sizeof(uchar2)*input_size));
	checkCudaErrors(cudaMalloc((void **) &d_output, sizeof(float2)*output_size));
	checkCudaErrors(cudaMalloc((void **) &d_coeff,  sizeof(float)*coeff_size));
	timer.Stop();
	if (DEBUG) printf("done in %g ms.", timer.Elapsed());

	//---------> Transfer coefficients to the device
	if (DEBUG) printf("\nCopy coefficients from host to device...\t");
	timer.Start();
	checkCudaErrors(cudaMemcpy(d_coeff, coeff, coeff_size*sizeof(float), cudaMemcpyHostToDevice));

	timer.Stop();
	transfer_in+=timer.Elapsed();
	if (DEBUG) printf("done in %g ms.\n", timer.Elapsed());
	
	//---------> Polyphase filter
	for (int r = 0; r < nRepeats; r++){
		//-----> Copy chunk of input data to a device
		timer.Start();
		checkCudaErrors(cudaMemcpy(d_input, &input[r*output_size], input_size*sizeof(uchar2), cudaMemcpyHostToDevice));
		timer.Stop();
		transfer_in+=timer.Elapsed();
		
		//-----> Compute Polyphase on the chunk
		Polyphase_GPU_init();
		Polyphase_GPU_benchmark(d_input, d_output, d_coeff, nChannels, nTaps, maxColumns, &fir_time, &fft_time);
		
		//-----> Copy chunk of output data to host
		timer.Start();
		checkCudaErrors(cudaMemcpy(&output[r*output_size], d_output, output_size*sizeof(float2), cudaMemcpyDeviceToHost));
		timer.Stop();
		transfer_out+=timer.Elapsed();
	}
	if (Sremainder>0){
		//-----> Copy chunk of input data to a device
		timer.Start();
		checkCudaErrors(cudaMemcpy(d_input, &input[nRepeats*output_size], (Sremainder+TAPS-1)*nChannels*sizeof(uchar2), cudaMemcpyHostToDevice));
		timer.Stop();
		transfer_in+=timer.Elapsed();
	
		//-----> Compute Polyphase on the chunk
		Polyphase_GPU_init();
		Polyphase_GPU_benchmark(d_input, d_output, d_coeff, nChannels, nTaps, Spectra_to_run, &fir_time, &fft_time);
		
		//-----> Copy chunk of output data to host
		timer.Start();
		checkCudaErrors(cudaMemcpy( &output[nRepeats*output_size], d_output, Sremainder*nChannels*sizeof(float2), cudaMemcpyDeviceToHost));
		timer.Stop();
		transfer_out+=timer.Elapsed();
	}


	//---------> error check -----
	checkCudaErrors(cudaGetLastError());
	
	//---------> Feeing allocated resources
	checkCudaErrors(cudaFree(d_input));
	checkCudaErrors(cudaFree(d_coeff));
	checkCudaErrors(cudaFree(d_output));


	if (DEBUG) printf("Number of spectra: %d;\nNumber of Channels: %d;\nNumber of Taps: %d;\nFIR filter execution time: %0.3f ms;\ncuFFT execution time: %0.3f ms;\nPolyphase execution time: %0.3f ms;\nData transfer time %0.3f ms\n",nSpectra,nChannels,nTaps, fir_time, fft_time, fir_time + fft_time, transfer_in + transfer_out);
	
	if (WRITE){ 
		char str[200];
		sprintf(str,"GPU-polyphase.dat");
		if (DEBUG) printf("\n Write results into file...\t");
		save_time(str, nSpectra, fir_time, fft_time, transfer_in, transfer_out, nChannels, nTaps, 0, THREADS_PER_BLOCK, DATA_SIZE);
		if (DEBUG) printf("\t done.\n-------------------------------------\n");
	}
	
}
