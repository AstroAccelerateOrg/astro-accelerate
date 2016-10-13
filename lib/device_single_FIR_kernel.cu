// Added by Karel Adamek 

#ifndef SINGLE_FIR_KERNEL_H_
#define SINGLE_FIR_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "AstroAccelerate/params.h"

<<<<<<< HEAD
__global__ void PD_FIR_GPU(float const* __restrict__ d_input, float *d_output, int nTaps, int nLoops, int nTimesamples) {
=======
__global__ void PD_FIR_GPU(float const* __restrict__ d_input, float *d_output, int nTaps, int nLoops, int nTimesamples)
{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	extern __shared__ float s_input[];

	int itemp,pos;
	float sum[PD_FIR_NWINDOWS];

<<<<<<< HEAD

	//----------------------------------------------
	//---- Reading data
	itemp=PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + nTaps - 1;
	for(int i=0; i<nLoops; i++){
		pos=i*PD_FIR_ACTIVE_WARPS*WARP + threadIdx.x;
		if(pos<itemp){
			s_input[pos]=d_input[blockIdx.y*nTimesamples + blockIdx.x*PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + pos];
		}
	}
	
	__syncthreads();
	
	
=======
	//----------------------------------------------
	//---- Reading data
	itemp = PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + nTaps - 1;
	for(int i=0; i<nLoops; i++)
	{
		pos = i*PD_FIR_ACTIVE_WARPS*WARP + threadIdx.x;
		if(pos<itemp)
			s_input[pos]=d_input[blockIdx.y*nTimesamples + blockIdx.x*PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + pos];
	}
	
	__syncthreads();
		
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	//----------------------------------------------
	//---- Calculating FIR version 2
	
	pos=PD_FIR_NWINDOWS*threadIdx.x;
<<<<<<< HEAD
	sum[0]=0;
	for(int t=0; t<nTaps; t++){
		sum[0]+=s_input[pos + t];
	}
	for(int i=1; i<PD_FIR_NWINDOWS; i++){
		pos=PD_FIR_NWINDOWS*threadIdx.x + i - 1;
		sum[i]=sum[i-1] - s_input[pos]+s_input[pos + nTaps];
=======
	sum[0] = 0;
	for(int t=0; t<nTaps; t++)
		sum[0] += s_input[pos + t];
	for(int i=1; i<PD_FIR_NWINDOWS; i++)
	{
		pos = PD_FIR_NWINDOWS*threadIdx.x + i - 1;
		sum[i] = sum[i-1] - s_input[pos]+s_input[pos + nTaps];
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	}

	//----------------------------------------------
	//---- Writing data	
<<<<<<< HEAD
	for(int i=0; i<PD_FIR_NWINDOWS; i++){
		pos=PD_FIR_NWINDOWS*threadIdx.x + i;
		d_output[blockIdx.y*nTimesamples + blockIdx.x*PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + pos]=sum[i];
	}
}


=======
	for(int i=0; i<PD_FIR_NWINDOWS; i++)
	{
		pos = PD_FIR_NWINDOWS*threadIdx.x + i;
		d_output[blockIdx.y*nTimesamples + blockIdx.x*PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + pos] = sum[i];
	}
}

>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
#endif