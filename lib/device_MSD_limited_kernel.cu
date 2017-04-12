// Added by Karel Adamek 

#ifndef MSD_LIMITED_KERNEL_H_
#define MSD_LIMITED_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

__global__ void MSD_GPU_limited(float const* __restrict__ d_input, float *d_output, int x_steps, int nColumns, int offset) {
	extern __shared__ float par[];

	int warp_id, local_id, dim_y, pos_x, pos_y, pos;
	float x;
	float M;
	float S;
	float j, jv;
	float ftemp;

	local_id = threadIdx.x & ( WARP - 1 );
	warp_id = threadIdx.x >> 5;
	dim_y = blockDim.x >> 5;

	//                           y              +                      x
	pos_y = ( blockIdx.y*dim_y + warp_id )*nColumns;
	pos_x = blockIdx.x*WARP*x_steps + local_id;
	// Since I assume that nRest >64 then I do not need any branching here.
	M = __ldg(&d_input[pos_y + pos_x]);
	S = 0;
	j = 1.0f;
	for (int xf = 1; xf < x_steps; xf++) {
		pos_x = pos_x + WARP;
		if (pos_x < ( nColumns - offset )) {
			x = __ldg(&d_input[pos_y + pos_x]);
			j = j + 1.0f;
			M = M + x;
			ftemp = ( j*x - M );
			S = S + 1.0f / ( j*( j - 1.0f ) )*ftemp*ftemp;
		}
	}

	par[threadIdx.x] = M;
	par[blockDim.x + threadIdx.x] = S;
	par[2*blockDim.x + threadIdx.x] = j;

	__syncthreads();

	// now all threads had saved their work, reduction follows

	for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1) {
		if (threadIdx.x < i) {
			// again since nRest>64 I do not need to check if j>0
			jv = par[2*blockDim.x + i + threadIdx.x];
			ftemp = ( jv / j*M - par[i + threadIdx.x] );
			S = S + par[blockDim.x + i + threadIdx.x] + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
			M = M + par[i + threadIdx.x];
			j = j + jv;

			par[threadIdx.x] = M;
			par[blockDim.x + threadIdx.x] = S;
			par[2*blockDim.x + threadIdx.x] = j;
		}

		__syncthreads();
	}

	// by now we should have only 32 partial results. shuffle reduction follows
	for (int q = HALF_WARP; q > 0; q = q >> 1) {
		jv = __shfl_down(j, q);
		ftemp = ( jv / j*M - __shfl_down(M, q) );
		S = S + __shfl_down(S, q) + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
		M = M + __shfl_down(M, q);
		j = j + jv;
	}

	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		pos = blockIdx.y*gridDim.x + blockIdx.x;
		d_output[3*pos] = M;
		d_output[3*pos + 1] = S;
		d_output[3*pos + 2] = j;
	}
}

__global__ void MSD_GPU_limited_final(float *d_input, float *d_output, int size) {
	__shared__ float Ms[WARP*WARP];
	__shared__ float Ss[WARP*WARP];
	__shared__ float js[WARP*WARP]; // I'll leave this as it is since this kernel launches in one copy thus I do not care how much shared memory it will eat up.

	// int warp_id;
	int pos;
	float M;
	float S;
	float j, jv;
	float ftemp;

	//warp_id = threadIdx.x >> 5;

	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos = threadIdx.x;
	if (size > blockDim.x) {
		M = d_input[3*pos];
		S = d_input[3*pos + 1];
		j = d_input[3*pos + 2];
		
		pos = pos + blockDim.x;
		while (pos < size) {
			jv = d_input[3*pos + 2];
			ftemp = ( jv / j*M - d_input[3*pos] );
			S = S + d_input[3*pos + 1] + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
			M = M + d_input[3*pos];
			j = j + jv;
			pos = pos + blockDim.x;
		}

		__syncthreads();

		Ms[threadIdx.x] = M;
		Ss[threadIdx.x] = S;
		js[threadIdx.x] = j;

		// now all threads had saved their work, reduction follows

		// first we must load initial values
		for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1) {
			if (threadIdx.x < i) {
				jv = js[i + threadIdx.x];
				ftemp = ( jv / j*M - Ms[i + threadIdx.x] );
				S = S + Ss[i + threadIdx.x] + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
				M = M + Ms[i + threadIdx.x];
				j = j + jv;

				Ms[threadIdx.x] = M;
				Ss[threadIdx.x] = S;
				js[threadIdx.x] = j;
			}
			__syncthreads();
		}

		// by now we should have only 32 partial results. shuffle reduction follows
		for (int q = HALF_WARP; q > 0; q = q >> 1)
		{
			jv = __shfl_down(j, q);
			ftemp = ( jv / j*M - __shfl_down(M, q) );
			S = S + __shfl_down(S, q) + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
			M = M + __shfl_down(M, q);
			j = j + jv;
		}

	}
	else {
		if (threadIdx.x == 0) {
			pos = 0;
			M = d_input[3*pos];
			S = d_input[3*pos + 1];
			j = d_input[3*pos + 2];
			
			for (pos = 1; pos < size; pos++) {
				jv = d_input[3*pos + 2];
				ftemp = ( jv / j*M - d_input[3*pos] );
				S = S + d_input[3*pos + 1] + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
				M = M + d_input[3*pos];
				j = j + jv;
			}
		}
	}

	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		d_output[0] = M / j;
		d_output[1] = sqrt(S / j);
		d_output[2] = j;

	}
}

__global__ void MSD_GPU_LA_ALL(float const* __restrict__ d_input, float *d_output, float *d_output_taps, int y_steps, int nTaps, int nTimesamples, int offset) {
	// This kernel will calculate only sigma_T I might merge it with sigma_1 if performing well
	// Shared memory required
	__shared__ float s_input[3*PD_NTHREADS];
	__shared__ float s_base[3*PD_NTHREADS];
	
	// MSD variables
	float M, S, j, jv;
	float M_b, S_b, j_b;
	float ftemp;
	// FIR variables
	int d, spos, local_id;
	size_t gpos;
	float2 ftemp1, ftemp2, ftemp3;
	float Bw[2];

	// First y coordinate is separated
	//-------------------> FIR
	spos = blockIdx.x*(2*PD_NTHREADS-nTaps) + 2*threadIdx.x;
	gpos = blockIdx.y*y_steps*nTimesamples + spos;
	Bw[0]=0; Bw[1]=0; j=0; j_b=0;
	if( (spos+4)<(nTimesamples-offset) ){
		ftemp1.x=d_input[gpos];	
		ftemp1.y=d_input[gpos+1];
		ftemp2.x=d_input[gpos+2];
		ftemp2.y=d_input[gpos+3];
		ftemp3.x=d_input[gpos+4];

		Bw[0]=ftemp1.x + ftemp1.y + ftemp2.x + ftemp2.y;
		Bw[1]=ftemp1.y + ftemp2.x + ftemp2.y + ftemp3.x;
		
		// Initialization of MSD variables
		M_b = ftemp1.x;
		S_b = 0;
		j_b = 1.0f;
		
		j_b = j_b + 1.0f;
		M_b = M_b + ftemp1.y;
		ftemp = ( j_b*ftemp1.y - M_b );
		S_b = S_b + 1.0f / ( j_b*( j_b - 1.0f ) )*ftemp*ftemp;
	}
	
	s_input[2*threadIdx.x] = Bw[0];
	s_input[2*threadIdx.x+1] = Bw[1];
	
	for(d=4; d<nTaps; d=d+4){
		__syncthreads();
		local_id = threadIdx.x-(d>>1);
		if(local_id>=0){
			Bw[0] = Bw[0] + s_input[2*local_id]; Bw[1] = Bw[1] + s_input[2*local_id+1];
		}
	}
	//-------------------> FIR
	
	// Initialization of MSD variables
	M = Bw[0];
	S = 0;
	j = 1.0f;
	// First addition (second actually, but first done this way) 
	j = j + 1.0f;
	M = M + Bw[1];
	ftemp = ( j*Bw[1] - M );
	S = S + 1.0f / ( j*( j - 1.0f ) )*ftemp*ftemp;
	
	// Rest of the iteration in y direction	
	for (int yf = 1; yf < y_steps; yf++) {
		__syncthreads();
		//-------------------> FIR
		spos = blockIdx.x*(2*PD_NTHREADS-nTaps) + 2*threadIdx.x;
		gpos = blockIdx.y*y_steps*nTimesamples + yf*nTimesamples + spos;
		Bw[0]=0; Bw[1]=0;
		if( (spos+4)<(nTimesamples-offset) ){
			ftemp1.x=d_input[gpos];	
			ftemp1.y=d_input[gpos+1];
			ftemp2.x=d_input[gpos+2];
			ftemp2.y=d_input[gpos+3];
			ftemp3.x=d_input[gpos+4];

			Bw[0]=ftemp1.x + ftemp1.y + ftemp2.x + ftemp2.y;
			Bw[1]=ftemp1.y + ftemp2.x + ftemp2.y + ftemp3.x;
			
			j_b = j_b + 1.0f;
			M_b = M_b + ftemp1.x;
			ftemp = ( j_b*ftemp1.x - M_b );
			S_b = S_b + 1.0f / ( j_b*( j_b - 1.0f ) )*ftemp*ftemp;
			
			j_b = j_b + 1.0f;
			M_b = M_b + ftemp1.y;
			ftemp = ( j_b*ftemp1.y - M_b );
			S_b = S_b + 1.0f / ( j_b*( j_b - 1.0f ) )*ftemp*ftemp;
		}
		
		s_input[2*threadIdx.x] = Bw[0];
		s_input[2*threadIdx.x+1] = Bw[1];
	
		for(d=4; d<nTaps; d=d+4){
			__syncthreads();
			local_id = threadIdx.x-(d>>1);
			if(local_id>=0){
				Bw[0] = Bw[0] + s_input[2*local_id]; Bw[1] = Bw[1] + s_input[2*local_id+1];
			}
		}
		//-------------------> FIR
		j = j + 1.0f;
		M = M + Bw[0];
		ftemp = ( j*Bw[0] - M );
		S = S + 1.0f / ( j*( j - 1.0f ) )*ftemp*ftemp;
	
		j = j + 1.0f;
		M = M + Bw[1];
		ftemp = ( j*Bw[1] - M );
		S = S + 1.0f / ( j*( j - 1.0f ) )*ftemp*ftemp;
	}
	
	__syncthreads();
	
	s_input[threadIdx.x] = 0;
	s_input[blockDim.x + threadIdx.x] = 0;
	s_input[2*blockDim.x + threadIdx.x] = 0;
	
	s_base[threadIdx.x] = 0;
	s_base[blockDim.x + threadIdx.x] = 0;
	s_base[2*blockDim.x + threadIdx.x] = 0;
	
	__syncthreads();
	
	if(local_id>=0) {
		// correct results have threads with 0 <= local_id < (2*PD_NTHREADS-nTaps); All other threads needs to store j=0  
		spos=blockIdx.x*(2*PD_NTHREADS-nTaps) + 2*threadIdx.x;
		if( spos<(nTimesamples-offset-nTaps) ) {
			s_input[local_id] = M;
			s_input[blockDim.x + local_id] = S;
			s_input[2*blockDim.x + local_id] = j;
			
			s_base[local_id] = M_b;
			s_base[blockDim.x + local_id] = S_b;
			s_base[2*blockDim.x + local_id] = j_b;
		}

	}
	__syncthreads();
	
	M=s_input[threadIdx.x];
	S=s_input[blockDim.x + threadIdx.x];
	j=s_input[2*blockDim.x + threadIdx.x];	
	
	M_b=s_base[threadIdx.x];
	S_b=s_base[blockDim.x + threadIdx.x];
	j_b=s_base[2*blockDim.x + threadIdx.x];	
	
	// We have finished loading data and creating partials for MSD and now we can reduce the MSD variables
	
	//------------------------------------------------------------------------------------
	//---------> StrDev of processed input
	for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1) {
		if (threadIdx.x < i) {
			jv = s_input[2*blockDim.x + i + threadIdx.x];
			if( ((int) jv)!=0){
				if(j==0){
					S = s_input[blockDim.x + i + threadIdx.x];
					M = s_input[i + threadIdx.x];
					j = jv;
				}
				else {
					ftemp = ( jv / j*M - s_input[i + threadIdx.x] );
					S = S + s_input[blockDim.x + i + threadIdx.x] + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
					M = M + s_input[i + threadIdx.x];
					j = j + jv;
				}
			}
			
			s_input[threadIdx.x] = M;
			s_input[blockDim.x + threadIdx.x] = S;
			s_input[2*blockDim.x + threadIdx.x] = j;
		}
		__syncthreads();		

	}


	// by now we should have only 32 partial results. shuffle reduction follows
	for (int q = HALF_WARP; q > 0; q = q >> 1) {
		jv = __shfl_down(j, q);
		if(jv!=0){
			if(j==0) {
				S = __shfl_down(S, q);
				M = __shfl_down(M, q);
				j = jv;
			}
			else {
				ftemp = ( jv / j*M - __shfl_down(M, q) );
				S = S + __shfl_down(S, q) + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
				M = M + __shfl_down(M, q);
				j = j + jv;
			}
		}
	}
	
	
	//------------------------------------------------------------------------------------
	//---------> StrDev of unprocessed input
	for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1) {
		if (threadIdx.x < i) {
			jv = s_base[2*blockDim.x + i + threadIdx.x];
			if( ((int) jv)!=0){
				if(j_b==0){
					S_b = s_base[blockDim.x + i + threadIdx.x];
					M_b = s_base[i + threadIdx.x];
					j_b = jv;
				}
				else {
					ftemp = ( jv / j_b*M_b - s_base[i + threadIdx.x] );
					S_b = S_b + s_base[blockDim.x + i + threadIdx.x] + ( j_b / ( jv*( j_b + jv ) ) )*ftemp*ftemp;
					M_b = M_b + s_base[i + threadIdx.x];
					j_b = j_b + jv;
				}
			}
			
			s_base[threadIdx.x] = M_b;
			s_base[blockDim.x + threadIdx.x] = S_b;
			s_base[2*blockDim.x + threadIdx.x] = j_b;
		}
		__syncthreads();		

	}


	// by now we should have only 32 partial results. shuffle reduction follows
	for (int q = HALF_WARP; q > 0; q = q >> 1) {
		jv = __shfl_down(j_b, q);
		if(jv!=0){
			if(j_b==0) {
				S_b = __shfl_down(S_b, q);
				M_b = __shfl_down(M_b, q);
				j_b = jv;
			}
			else {
				ftemp = ( jv / j_b*M_b - __shfl_down(M_b, q) );
				S_b = S_b + __shfl_down(S_b, q) + ( j_b / ( jv*( j_b + jv ) ) )*ftemp*ftemp;
				M_b = M_b + __shfl_down(M_b, q);
				j_b = j_b + jv;
			}
		}
	}

	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		gpos = blockIdx.y*gridDim.x + blockIdx.x;
		d_output_taps[3*gpos] = M;
		d_output_taps[3*gpos + 1] = S;
		d_output_taps[3*gpos + 2] = j;
		
		d_output[3*gpos] = M_b;
		d_output[3*gpos + 1] = S_b;
		d_output[3*gpos + 2] = j_b;
	}
}


__global__ void MSD_GPU_limited_final_create_linear_approx(float *d_input, float *d_output, float *d_MSD_base, int nTaps, int size) {
	__shared__ float Ms[WARP*WARP];
	__shared__ float Ss[WARP*WARP];
	__shared__ float js[WARP*WARP]; // I'll leave this as it is since this kernel launches in one copy thus I do not care how much shared memory it will eat up.

	// int warp_id;
	int pos;
	float M;
	float S;
	float j, jv;
	float ftemp;

	//warp_id = threadIdx.x >> 5;

	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos = threadIdx.x;
	if (size > blockDim.x) {
		M = d_input[3*pos];
		S = d_input[3*pos + 1];
		j = d_input[3*pos + 2];
		
		pos = pos + blockDim.x;
		while (pos < size) {
			jv = d_input[3*pos + 2];
			ftemp = ( jv / j*M - d_input[3*pos] );
			S = S + d_input[3*pos + 1] + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
			M = M + d_input[3*pos];
			j = j + jv;
			pos = pos + blockDim.x;
		}

		__syncthreads();

		Ms[threadIdx.x] = M;
		Ss[threadIdx.x] = S;
		js[threadIdx.x] = j;

		// now all threads had saved their work, reduction follows

		// first we must load initial values
		for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1) {
			if (threadIdx.x < i) {
				jv = js[i + threadIdx.x];
				ftemp = ( jv / j*M - Ms[i + threadIdx.x] );
				S = S + Ss[i + threadIdx.x] + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
				M = M + Ms[i + threadIdx.x];
				j = j + jv;

				Ms[threadIdx.x] = M;
				Ss[threadIdx.x] = S;
				js[threadIdx.x] = j;
			}
			__syncthreads();
		}

		// by now we should have only 32 partial results. shuffle reduction follows
		for (int q = HALF_WARP; q > 0; q = q >> 1)
		{
			jv = __shfl_down(j, q);
			ftemp = ( jv / j*M - __shfl_down(M, q) );
			S = S + __shfl_down(S, q) + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
			M = M + __shfl_down(M, q);
			j = j + jv;
		}

	}
	else {
		if (threadIdx.x == 0) {
			pos = 0;
			M = d_input[3*pos];
			S = d_input[3*pos + 1];
			j = d_input[3*pos + 2];
			;
			for (pos = 1; pos < size; pos++) {
				jv = d_input[3*pos + 2];
				ftemp = ( jv / j*M - d_input[3*pos] );
				S = S + d_input[3*pos + 1] + ( j / ( jv*( j + jv ) ) )*ftemp*ftemp;
				M = M + d_input[3*pos];
				j = j + jv;
			}
		}
	}

	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		d_output[0] = d_MSD_base[0];//M / j;
		d_output[1] = d_MSD_base[1];//sqrt(S / j)
		d_output[2] = (sqrt(S / j) - d_MSD_base[1])/( (float) (nTaps-1));
	}
}

#endif
