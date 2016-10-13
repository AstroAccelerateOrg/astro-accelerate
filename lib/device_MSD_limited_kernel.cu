// Added by Karel Adamek 

#ifndef MSD_LIMITED_KERNEL_H_
#define MSD_LIMITED_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "AstroAccelerate/params.h"



<<<<<<< HEAD
__global__ void MSD_GPU_limited(float const* __restrict__ d_input, float *d_output, int x_steps, int nColumns, int offset) {
=======
__global__ void MSD_GPU_limited(float const* __restrict__ d_input, float *d_output, int x_steps, int nColumns, int offset)
{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	extern __shared__ float par[];
	
	int warp_id, local_id, dim_y, pos_x, pos_y, pos;
	float x;
	float M;
	float S;
	float j,jv;
	float ftemp;
	
	local_id = threadIdx.x & (WARP - 1);
	warp_id = threadIdx.x>>5;
	dim_y = blockDim.x>>5;
	
	//                           y              +                      x
	pos_y = (blockIdx.y*dim_y + warp_id)*nColumns;
	pos_x = blockIdx.x*WARP*x_steps + local_id;
	// Since I assume that nRest >64 then I do not need any branching here.
<<<<<<< HEAD
	M=__ldg(&d_input[pos_y + pos_x]);
	S=0;
	j=1.0f;
	for(int xf=1; xf<x_steps; xf++){
		pos_x = pos_x + WARP;
		if(pos_x<(nColumns-offset)){
=======
	M = __ldg(&d_input[pos_y + pos_x]);
	S = 0;
	j = 1.0f;
	for(int xf=1; xf<x_steps; xf++)
	{
		pos_x = pos_x + WARP;
		if(pos_x<(nColumns-offset))
		{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
			x = __ldg(&d_input[pos_y + pos_x]);
			j = j+1.0f;
			M = M + x;
			ftemp = (j*x - M);
			S = S + 1.0f/(j*(j-1.0f))*ftemp*ftemp;
		}
	}
	
<<<<<<< HEAD
	par[threadIdx.x]=M;
	par[blockDim.x + threadIdx.x]=S;
	par[2*blockDim.x + threadIdx.x]=j;
=======
	par[threadIdx.x] = M;
	par[blockDim.x + threadIdx.x] = S;
	par[2*blockDim.x + threadIdx.x] = j;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	
	__syncthreads();
	
	// now all threads had saved their work, reduction follows
	
<<<<<<< HEAD
	for(int i=(blockDim.x>>1); i>HALF_WARP; i=i>>1){
		if(threadIdx.x<i){
			// again since nRest>64 I do not need to check if j>0
			jv=par[2*blockDim.x + i + threadIdx.x];
			ftemp = (jv/j*M - par[i + threadIdx.x]);
			S = S + par[blockDim.x + i + threadIdx.x] + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + par[i + threadIdx.x];
			j=j+jv;
			
			par[threadIdx.x]=M;
			par[blockDim.x + threadIdx.x]=S;
			par[2*blockDim.x + threadIdx.x]=j;
		}
		
=======
	for(int i=(blockDim.x>>1); i>HALF_WARP; i=i>>1)
	{
		if(threadIdx.x<i)
		{
			// again since nRest>64 I do not need to check if j>0
			jv = par[2*blockDim.x + i + threadIdx.x];
			ftemp = (jv/j*M - par[i + threadIdx.x]);
			S = S + par[blockDim.x + i + threadIdx.x] + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + par[i + threadIdx.x];
			j = j+jv;
			
			par[threadIdx.x] = M;
			par[blockDim.x + threadIdx.x] = S;
			par[2*blockDim.x + threadIdx.x] = j;
		}		
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		__syncthreads();
	}
	
	
	// by now we should have only 32 partial results. shuffle reduction follows
<<<<<<< HEAD
	for(int q=HALF_WARP; q>0; q=q>>1){
		jv=__shfl_down(j, q);
		ftemp = (jv/j*M - __shfl_down(M, q));
		S = S + __shfl_down(S, q) + (j/(jv*(j+jv)))*ftemp*ftemp;
		M = M + __shfl_down(M, q);
		j=j+jv;
=======
	for(int q=HALF_WARP; q>0; q=q>>1)
	{
		jv = __shfl_down(j, q);
		ftemp = (jv/j*M - __shfl_down(M, q));
		S = S + __shfl_down(S, q) + (j/(jv*(j+jv)))*ftemp*ftemp;
		M = M + __shfl_down(M, q);
		j = j+jv;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	}
	
	//----------------------------------------------
	//---- Writing data
<<<<<<< HEAD
	if(threadIdx.x==0){
=======
	if(threadIdx.x == 0)
	{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		pos = blockIdx.y*gridDim.x + blockIdx.x;
		d_output[3*pos] = M;
		d_output[3*pos + 1] = S;
		d_output[3*pos + 2] = j;
	}
}

<<<<<<< HEAD


__global__ void MSD_GPU_limited_final(float *d_input, float *d_output, int size) {
=======
__global__ void MSD_GPU_limited_final(float *d_input, float *d_output, int size)
{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	__shared__ float Ms[WARP*WARP];
	__shared__ float Ss[WARP*WARP];
	__shared__ float js[WARP*WARP]; // I'll leave this as it is since this kernel launches in one copy thus I do not care how much shared memory it will eat up.
	
	int warp_id, pos;
	float M;
	float S;
	float j, jv;
	float ftemp;
	
	warp_id = threadIdx.x>>5;
	
<<<<<<< HEAD
	
	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos=threadIdx.x;
	if(size>blockDim.x){
		M=d_input[3*pos];
		S=d_input[3*pos+1];
		j=d_input[3*pos+2];;
		pos = pos + blockDim.x;
		while (pos<size){
			jv=d_input[3*pos+2];
			ftemp = (jv/j*M - d_input[3*pos]);
			S = S + d_input[3*pos+1] + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + d_input[3*pos];
			j=j+jv;
=======
	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos = threadIdx.x;
	if(size>blockDim.x)
	{
		M = d_input[3*pos];
		S = d_input[3*pos+1];
		j = d_input[3*pos+2];;
		pos = pos + blockDim.x;
		while (pos<size)
		{
			jv = d_input[3*pos+2];
			ftemp = (jv/j*M - d_input[3*pos]);
			S = S + d_input[3*pos+1] + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + d_input[3*pos];
			j = j+jv;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
			pos = pos + blockDim.x;
		}
		
		__syncthreads();
		
<<<<<<< HEAD
		Ms[threadIdx.x]=M;
		Ss[threadIdx.x]=S;
		js[threadIdx.x]=j;
=======
		Ms[threadIdx.x] = M;
		Ss[threadIdx.x] = S;
		js[threadIdx.x] = j;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		
		// now all threads had saved their work, reduction follows
		
		// first we must load initial values
<<<<<<< HEAD
		for(int i=(blockDim.x>>1); i>HALF_WARP; i=i>>1){
			if(threadIdx.x<i){
				jv=js[i + threadIdx.x];
				ftemp = (jv/j*M - Ms[i + threadIdx.x]);
				S = S + Ss[i + threadIdx.x] + (j/(jv*(j+jv)))*ftemp*ftemp;
				M = M + Ms[i + threadIdx.x];
				j=j+jv;
				
				Ms[threadIdx.x]=M;
				Ss[threadIdx.x]=S;
				js[threadIdx.x]=j;
=======
		for(int i=(blockDim.x>>1); i>HALF_WARP; i=i>>1)
		{
			if(threadIdx.x<i)
			{
				jv = js[i + threadIdx.x];
				ftemp = (jv/j*M - Ms[i + threadIdx.x]);
				S = S + Ss[i + threadIdx.x] + (j/(jv*(j+jv)))*ftemp*ftemp;
				M = M + Ms[i + threadIdx.x];
				j = j+jv;
				
				Ms[threadIdx.x] = M;
				Ss[threadIdx.x] = S;
				js[threadIdx.x] = j;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
			}
			__syncthreads();
		}
		
		// by now we should have only 32 partial results. shuffle reduction follows
<<<<<<< HEAD
		for(int q=HALF_WARP; q>0; q=q>>1){
			jv=__shfl_down(j, q);
			ftemp = (jv/j*M - __shfl_down(M, q));
			S = S + __shfl_down(S, q) + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + __shfl_down(M, q);
			j=j+jv;
=======
		for(int q=HALF_WARP; q>0; q=q>>1)
		{
			jv = __shfl_down(j, q);
			ftemp = (jv/j*M - __shfl_down(M, q));
			S = S + __shfl_down(S, q) + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + __shfl_down(M, q);
			j = j+jv;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		}
		
	}
	else {
<<<<<<< HEAD
		if(threadIdx.x==0){
			pos=0;
			M=d_input[3*pos];
			S=d_input[3*pos+1];
			j=d_input[3*pos+2];;
			for(pos=1; pos<size; pos++){
				jv=d_input[3*pos+2];
				ftemp = (jv/j*M - d_input[3*pos]);
				S = S + d_input[3*pos+1] + (j/(jv*(j+jv)))*ftemp*ftemp;
				M = M + d_input[3*pos];
				j=j+jv;
=======
		if(threadIdx.x==0)
		{
			pos = 0;
			M = d_input[3*pos];
			S = d_input[3*pos+1];
			j = d_input[3*pos+2];;
			for(pos=1; pos<size; pos++)
			{
				jv = d_input[3*pos+2];
				ftemp = (jv/j*M - d_input[3*pos]);
				S = S + d_input[3*pos+1] + (j/(jv*(j+jv)))*ftemp*ftemp;
				M = M + d_input[3*pos];
				j = j+jv;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
			}
		}
	}
	
	//----------------------------------------------
	//---- Writing data
<<<<<<< HEAD
	if(threadIdx.x==0){
		d_output[0] = M/j;
		d_output[1] = sqrt(S/j);
		d_output[2] = j;
		
	}
}











=======
	if(threadIdx.x==0)
	{
		d_output[0] = M/j;
		d_output[1] = sqrt(S/j);
		d_output[2] = j;		
	}
}

>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
#endif