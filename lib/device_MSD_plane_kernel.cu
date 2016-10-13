// Added by Karel Adamek 

#ifndef MSD_PLANE_KERNEL_H_
#define MSD_PLANE_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "AstroAccelerate/params.h"


<<<<<<< HEAD
__global__ void MSD_GPU(float2 const* __restrict__ d_input, float *d_output) {
=======
__global__ void MSD_GPU(float2 const* __restrict__ d_input, float *d_output)
{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	__shared__ float Ms[WARP*MSD_WARPS_PER_BLOCK];
	__shared__ float Ss[WARP*MSD_WARPS_PER_BLOCK];
	
	int warp_id, local_id, pos;
	float2 x;
	float M;
	float S;
	float j;
	float ftemp;
	
	local_id = threadIdx.x & (WARP - 1);
	warp_id = threadIdx.x>>5;
	
	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos=blockIdx.x*MSD_WARPS_PER_BLOCK*WARP*MSD_ELEM_PER_THREAD + warp_id*WARP*MSD_ELEM_PER_THREAD + local_id;
	x = __ldg(&d_input[pos]);
	M = x.x;
	S = 0;
	j = 1.0f;
	
	j = j+1.0f;
	M = M + x.y;
	ftemp = (j*x.y - M);
	S = S + 1.0f/(j*(j-1.0f))*ftemp*ftemp;
	
<<<<<<< HEAD
	for(int i=1; i<MSD_ELEM_PER_THREAD; i++){
=======
	for(int i=1; i<MSD_ELEM_PER_THREAD; i++)
	{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		pos = pos + WARP;
		x = __ldg(&d_input[pos]);
		
		j = j+1.0f;
		M = M + x.x;
		ftemp = (j*x.x - M);
		S = S + 1.0f/(j*(j-1.0f))*ftemp*ftemp;
		
		j = j+1.0f;
		M = M + x.y;
		ftemp = (j*x.y - M);
		S = S + 1.0f/(j*(j-1.0f))*ftemp*ftemp;
	}
<<<<<<< HEAD
	Ms[threadIdx.x]=M;
	Ss[threadIdx.x]=S;
=======
	Ms[threadIdx.x] = M;
	Ss[threadIdx.x] = S;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	
	__syncthreads();
	
	// now all threads had saved their work, reduction follows
	
	// first we must load initial values
	//j=2*MSD_ELEM_PER_THREAD; // value of j is preserved during kernel's execution
<<<<<<< HEAD
	for(int i=(blockDim.x>>1); i>HALF_WARP; i=i>>1){
		if(threadIdx.x<i){
			j=j*2;
=======
	for(int i=(blockDim.x>>1); i>HALF_WARP; i=i>>1)
	{
		if(threadIdx.x<i){
			j = j*2;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
			ftemp = (M - Ms[i + threadIdx.x]);
			S = S + Ss[i + threadIdx.x] + (1.0f/j)*ftemp*ftemp;
			M = M + Ms[i + threadIdx.x];
			
<<<<<<< HEAD
			Ms[threadIdx.x]=M;
			Ss[threadIdx.x]=S;
=======
			Ms[threadIdx.x] = M;
			Ss[threadIdx.x] = S;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		}
		__syncthreads();
	}
	
	// by now we should have only 32 partial results. shuffle reduction follows
<<<<<<< HEAD
	for(int q=HALF_WARP; q>0; q=q>>1){
		j=j*2;
=======
	for(int q=HALF_WARP; q>0; q=q>>1)
	{
		j = j*2;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		ftemp = (M - __shfl_down(M, q));
		S = S + __shfl_down(S, q) + (1.0f/j)*ftemp*ftemp;
		M = M + __shfl_down(M, q);
	}
	
	//----------------------------------------------
	//---- Writing data
<<<<<<< HEAD
	if(warp_id==0 && threadIdx.x==0){
=======
	if(warp_id == 0 && threadIdx.x == 0)
	{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		d_output[2*blockIdx.x] = M;
		d_output[2*blockIdx.x + 1] = S;
	}
}

<<<<<<< HEAD



__global__ void MSD_GPU_remainder(float const* __restrict__ d_input, float *d_output, int remainder) {
=======
__global__ void MSD_GPU_remainder(float const* __restrict__ d_input, float *d_output, int remainder)
{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	__shared__ float Ms[WARP*MSD_WARPS_PER_BLOCK];
	__shared__ float Ss[WARP*MSD_WARPS_PER_BLOCK];
	__shared__ float js[WARP*MSD_WARPS_PER_BLOCK];
	
	int warp_id, pos;
	float x;
	float M;
	float S;
	float j, jv;
	float ftemp;
	
	warp_id = threadIdx.x>>5;
	
	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos=threadIdx.x;
<<<<<<< HEAD
	if(remainder>blockDim.x){
		M=__ldg(&d_input[pos]);
		S=0;
		j=1.0f;
		pos = pos + blockDim.x;
		while (pos<remainder){
=======
	if(remainder>blockDim.x)
	{
		M = __ldg(&d_input[pos]);
		S = 0;
		j = 1.0f;
		pos = pos + blockDim.x;
		while (pos<remainder)
		{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
			x = __ldg(&d_input[pos]);
			j = j+1.0f;
			M = M + x;
			ftemp = (j*x - M);
			S = S + 1.0f/(j*(j-1.0f))*ftemp*ftemp;
			pos = pos + blockDim.x;
		}
		
<<<<<<< HEAD
		Ms[threadIdx.x]=M;
		Ss[threadIdx.x]=S;
		js[threadIdx.x]=j;
=======
		Ms[threadIdx.x] = M;
		Ss[threadIdx.x] = S;
		js[threadIdx.x] = j;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		
		__syncthreads();
		
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
			if(threadIdx.x<i){
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
		}
		
	}
	else {
		if(threadIdx.x==0){// This assumes remainder to be small < 32
			pos=0;
			M=__ldg(&d_input[pos]);
			S=0;
			j=1.0f;
			for(pos=1; pos<remainder; pos++){
=======
		for(int q=HALF_WARP; q>0; q=q>>1)
		{
			jv = __shfl_down(j, q);
			ftemp = (jv/j*M - __shfl_down(M, q));
			S = S + __shfl_down(S, q) + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + __shfl_down(M, q);
			j = j+jv;
		} 		
	}
	else 
	{
		if(threadIdx.x == 0)
		{// This assumes remainder to be small < 32
			pos = 0;
			M = __ldg(&d_input[pos]);
			S = 0;
			j = 1.0f;
			for(pos=1; pos<remainder; pos++)
			{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
				x = __ldg(&d_input[pos]);
				j = j+1.0f;
				M = M + x;
				ftemp = (j*x - M);
				S = S + 1.0f/(j*(j-1.0f))*ftemp*ftemp;
			}
		}
	}
	
	//----------------------------------------------
	//---- Writing data
<<<<<<< HEAD
	if(warp_id==0 && threadIdx.x==0){
=======
	if(warp_id == 0 && threadIdx.x == 0)
	{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		d_output[0] = M;
		d_output[1] = S;
		d_output[2] = j;
	}
}


<<<<<<< HEAD
__global__ void MSD_GPU_final(float *d_input, float *d_output, int size, int tail, float nElements) {
=======
__global__ void MSD_GPU_final(float *d_input, float *d_output, int size, int tail, float nElements)
{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	__shared__ float Ms[WARP*MSD_WARPS_PER_BLOCK];
	__shared__ float Ss[WARP*MSD_WARPS_PER_BLOCK];
	__shared__ float js[WARP*MSD_WARPS_PER_BLOCK];
	
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
		M=d_input[2*pos];
		S=d_input[2*pos+1];
		j=2*WARP*MSD_ELEM_PER_THREAD*MSD_WARPS_PER_BLOCK;
		jv=2*WARP*MSD_ELEM_PER_THREAD*MSD_WARPS_PER_BLOCK;
		pos = pos + blockDim.x;
		while (pos<size){
			ftemp = (jv/j*M - d_input[2*pos]);
			S = S + d_input[2*pos+1] + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + d_input[2*pos];
			j=j+jv;
=======
		
	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos=threadIdx.x;
	if(size>blockDim.x)
	{
		M = d_input[2*pos];
		S = d_input[2*pos+1];
		j = 2*WARP*MSD_ELEM_PER_THREAD*MSD_WARPS_PER_BLOCK;
		jv=  2*WARP*MSD_ELEM_PER_THREAD*MSD_WARPS_PER_BLOCK;
		pos = pos + blockDim.x;
		while (pos<size)
		{
			ftemp = (jv/j*M - d_input[2*pos]);
			S = S + d_input[2*pos+1] + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + d_input[2*pos];
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
		}
		
		if(tail>0 && threadIdx.x==0){
			jv=d_input[2*size+2];
=======
		for(int q=HALF_WARP; q>0; q=q>>1)
		{
			jv = __shfl_down(j, q);
			ftemp = (jv/j*M - __shfl_down(M, q));
			S = S + __shfl_down(S, q) + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + __shfl_down(M, q);
			j = j+jv;
		}
		
		if(tail>0 && threadIdx.x == 0)
		{
			jv = d_input[2*size+2];
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
			ftemp = (jv/j*M - d_input[2*size]);
			S = S + d_input[2*size+1] + (j/(jv*(j+jv)))*ftemp*ftemp;
			M = M + d_input[2*size];
			j=j+jv;
		}
	}
<<<<<<< HEAD
	else {
		if(threadIdx.x==0){
			pos=0;
			M=d_input[2*pos];
			S=d_input[2*pos+1];
			j=2*WARP*MSD_ELEM_PER_THREAD*MSD_WARPS_PER_BLOCK;
			for(pos=1; pos<size; pos++){
				j=j*2;
				ftemp = (M - d_input[2*pos]);
				S = S + d_input[2*pos+1] + (1.0f/j)*ftemp*ftemp;
				M = M + d_input[2*pos];
			}
			
			if(tail>0){
=======
	else 
	{
		if(threadIdx.x == 0)
		{
			pos = 0;
			M = d_input[2*pos];
			S = d_input[2*pos+1];
			j = 2*WARP*MSD_ELEM_PER_THREAD*MSD_WARPS_PER_BLOCK;
			for(pos=1; pos<size; pos++)
			{
				j = j*2;
				ftemp = (M - d_input[2*pos]);
				S = S + d_input[2*pos+1] + (1.0f/j)*ftemp*ftemp;
				M = M + d_input[2*pos];
			}			
			if(tail>0)
			{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
				jv=d_input[2*size+2];
				ftemp = (jv/j*M - d_input[2*size]);
				S = S + d_input[2*size+1] + (j/(jv*(j+jv)))*ftemp*ftemp;
				M = M + d_input[2*size];
<<<<<<< HEAD
				j=j+jv;
=======
				j = j+jv;
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
			}
		}
	}
	
	//----------------------------------------------
	//---- Writing data
<<<<<<< HEAD
	if(warp_id==0 && threadIdx.x==0){
		d_output[0] = M/nElements;
		d_output[1] = sqrt(S/nElements);
		d_output[2] = nElements;
		
	}
}



=======
	if(warp_id==0 && threadIdx.x == 0)
	{
		d_output[0] = M/nElements;
		d_output[1] = sqrt(S/nElements);
		d_output[2] = nElements;		
	}
}

>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
#endif