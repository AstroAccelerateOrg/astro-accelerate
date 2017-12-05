// Added by Karel Adamek 

#ifndef MSD_LIMITED_KERNEL_H_
#define MSD_LIMITED_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"



//----------------------------------------------------------------------------------------
//------------- Device functions
__device__ __inline__ void Initiate(float *M, float *S, float *j, float element){
	*M = element;
	*S = 0;
	*j = 1.0f;
}

__device__ __inline__ void Add_one(float *M, float *S, float *j, float element){
	float ftemp;
	*j = (*j) + 1.0f;
	*M = (*M) + element;
	ftemp = ( (*j)*element - (*M) );
	*S = (*S) + 1.0f / ( (*j)*( (*j) - 1.0f ) )*ftemp*ftemp;
}

__device__ __inline__ void Merge(float *A_M, float *A_S, float *A_j, float B_M, float B_S, float B_j){
	float ftemp;
	
	ftemp = ( B_j / (*A_j)*(*A_M) - B_M );
	(*A_S) = (*A_S) + B_S + ( (*A_j) / ( B_j*( (*A_j) + B_j ) ) )*ftemp*ftemp;
	(*A_M) = (*A_M) + B_M;
	(*A_j) = (*A_j) + B_j;
}

__device__ __inline__ void Reduce_SM(float *M, float *S, float *j, float *s_input){
	float jv;
	
	(*M)=s_input[threadIdx.x];
	(*S)=s_input[blockDim.x + threadIdx.x];
	(*j)=s_input[2*blockDim.x + threadIdx.x];
	
	for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1) {
		if (threadIdx.x < i) {
			jv = s_input[2*blockDim.x + i + threadIdx.x];
			if( ((int) jv)!=0){
				if( (*j)==0 ){
					(*S) = s_input[blockDim.x + i + threadIdx.x];
					(*M) = s_input[i + threadIdx.x];
					(*j) = jv;
				}
				else {
					Merge(M, S, j, s_input[i + threadIdx.x], s_input[blockDim.x + i + threadIdx.x], jv);
				}
			}
			
			s_input[threadIdx.x] = (*M);
			s_input[blockDim.x + threadIdx.x] = (*S);
			s_input[2*blockDim.x + threadIdx.x] = (*j);
		}
		__syncthreads();
	}
}

__device__ __inline__ void Reduce_SM_regular(float *M, float *S, float *j, float *s_input){
	(*M)=s_input[threadIdx.x];
	(*S)=s_input[blockDim.x + threadIdx.x];
	(*j)=s_input[2*blockDim.x + threadIdx.x];
	
	for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1) {
		if (threadIdx.x < i) {
			Merge(M, S, j, s_input[i + threadIdx.x], s_input[blockDim.x + i + threadIdx.x], s_input[2*blockDim.x + i + threadIdx.x]);
			
			s_input[threadIdx.x] = (*M);
			s_input[blockDim.x + threadIdx.x] = (*S);
			s_input[2*blockDim.x + threadIdx.x] = (*j);
		}
		__syncthreads();
	}
}

__device__ __inline__ void Reduce_WARP(float *M, float *S, float *j){
	float B_M, B_S, B_j;
	
	for (int q = HALF_WARP; q > 0; q = q >> 1) {
		B_M = __shfl_down((*M), q);
		B_S = __shfl_down((*S), q);
		B_j = __shfl_down((*j), q);
		
		if(B_j>0){
			if( (*j)==0 ) {
				(*S) = B_S;
				(*M) = B_M;
				(*j) = B_j;
			}
			else {
				Merge(M, S, j, B_M, B_S, B_j);
			}
		}
	}
}

__device__ __inline__ void Reduce_WARP_regular(float *M, float *S, float *j){
	for (int q = HALF_WARP; q > 0; q = q >> 1) {
		Merge(M, S, j, __shfl_down((*M), q), __shfl_down((*S), q), __shfl_down((*j), q));
	}
}

__device__ void Sum_partials_regular(float *M, float *S, float *j, float *d_input, float *s_input, int size){
	int pos;
	
	//----------------------------------------------
	//---- Summing partials
	pos = threadIdx.x;
	if (size > blockDim.x) {
		(*M) = d_input[3*pos];
		(*S) = d_input[3*pos + 1];
		(*j) = d_input[3*pos + 2];
		
		pos = pos + blockDim.x;
		while (pos < size) {
			Merge( M, S, j, d_input[3*pos], d_input[3*pos + 1], d_input[3*pos + 2]);
			pos = pos + blockDim.x;
		}

		s_input[threadIdx.x] = (*M);
		s_input[blockDim.x + threadIdx.x] = (*S);
		s_input[2*blockDim.x + threadIdx.x] = (*j);
		
		__syncthreads();

		Reduce_SM_regular( M, S, j, s_input);
		Reduce_WARP_regular(M, S, j);
	}
	else {
		if (threadIdx.x == 0) {
			pos = 0;
			(*M) = d_input[3*pos];
			(*S) = d_input[3*pos + 1];
			(*j) = d_input[3*pos + 2];
			
			for (pos = 1; pos < size; pos++) {
				Merge( M, S, j, d_input[3*pos], d_input[3*pos + 1], d_input[3*pos + 2]);
			}
		}
	}
	//---- Summing partials
	//----------------------------------------------
}

__device__ void Sum_partials_nonregular(float *M, float *S, float *j, float *d_input, float *s_input, int size){
	int pos;
	float jv;
	
	//----------------------------------------------
	//---- Summing partials
	pos = threadIdx.x;
	if (size > blockDim.x) {
		(*M) = 0;	(*S) = 0;	(*j) = 0;
		while (pos < size) {
			jv = d_input[3*pos + 2];
			if( ((int) jv)>0 ){
				if( (int) (*j)==0 ){
					(*M) = d_input[3*pos]; 
					(*S) = d_input[3*pos + 1];
					(*j) = jv;
				}
				else {
					Merge( M, S, j, d_input[3*pos], d_input[3*pos + 1], jv);
				}
			}
			pos = pos + blockDim.x;
		}

		s_input[threadIdx.x] = (*M);
		s_input[blockDim.x + threadIdx.x] = (*S);
		s_input[2*blockDim.x + threadIdx.x] = (*j);
		
		__syncthreads();

		Reduce_SM( M, S, j, s_input);
		Reduce_WARP(M, S, j);
	}
	else {
		if (threadIdx.x == 0) {
			pos = 0;
			(*M) = 0;	(*S) = 0;	(*j) = 0;
			for (pos = 1; pos < size; pos++) {
				jv = d_input[3*pos + 2];
				if( ((int) jv)!=0 ){
					if( (int) (*j)==0 ){
						(*M) = d_input[3*pos]; 
						(*S) = d_input[3*pos + 1];
						(*j) = jv;
					}
					else {
						Merge( M, S, j, d_input[3*pos], d_input[3*pos + 1], jv);
					}
				}
			}
		}
	}
	//---- Summing partials
	//----------------------------------------------
}

//------------- Device functions
//----------------------------------------------------------------------------------------


// Computes partials for mean and standard deviation of the data with offset at the end
// PD_THREADS could be replaced it is not required to be #defined
__global__ void MSD_GPU_limited(float const* __restrict__ d_input, float *d_output, int y_steps, int nTimesamples, int offset) {
	__shared__ float s_input[3*PD_NTHREADS];
	float M, S, j, ftemp;
	
	int spos = blockIdx.x*PD_NTHREADS + threadIdx.x;
	int gpos = blockIdx.y*y_steps*nTimesamples + spos;
	M=0;	S=0;	j=0;
	if( spos<(nTimesamples-offset) ){
		
		ftemp=__ldg(&d_input[gpos]);
		Initiate( &M, &S, &j, ftemp);
		
		gpos = gpos + nTimesamples;
		for (int yf = 1; yf < y_steps; yf++) {
			ftemp=__ldg(&d_input[gpos]);
			Add_one( &M, &S, &j, ftemp);
			gpos = gpos + nTimesamples;
		}
	}
	
	s_input[threadIdx.x] = M;
	s_input[blockDim.x + threadIdx.x] = S;
	s_input[2*blockDim.x + threadIdx.x] = j;
	
	__syncthreads();
	
	Reduce_SM( &M, &S, &j, s_input );
	Reduce_WARP( &M, &S, &j);
	
	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		gpos = blockIdx.y*gridDim.x + blockIdx.x;
		d_output[3*gpos] = M;
		d_output[3*gpos + 1] = S;
		d_output[3*gpos + 2] = j;
	}
}

// Computes mean and standard deviation from partial
__global__ void MSD_GPU_final_regular(float *d_input, float *d_output, int size) {
	__shared__ float s_input[3*WARP*WARP];

	float M, S, j;
	
	Sum_partials_regular( &M, &S, &j, d_input, s_input, size);

	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		d_output[0] = M / j;
		d_output[1] = sqrt(S / j);
		d_output[2] = j;
	}
}

__global__ void MSD_GPU_final_regular(float *d_input, float *d_MSD, float *d_pp, int size) {
	__shared__ float s_input[3*WARP*WARP];

	float M, S, j;
	
	Sum_partials_regular( &M, &S, &j, d_input, s_input, size);

	if(d_pp[2]>0){
		Merge(&M, &S, &j, d_pp[0], d_pp[1], d_pp[2]);
	}
	
	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		d_MSD[0] = M / j;
		d_MSD[1] = sqrt(S / j);
		d_MSD[2] = j;
		d_pp[0] = M;
		d_pp[1] = S;
		d_pp[2] = j;
	}
}


__global__ void MSD_GPU_final_nonregular(float *d_input, float *d_MSD, int size) {
	__shared__ float s_input[3*WARP*WARP];
	
	float M, S, j;

	Sum_partials_nonregular( &M, &S, &j, d_input, s_input, size);
	
	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		d_MSD[0] = M / j;
		d_MSD[1] = sqrt(S / j);
		d_MSD[2] = j;
	}
}

__global__ void MSD_GPU_final_nonregular(float *d_input, float *d_MSD, float *d_pp, int size) {
	__shared__ float s_input[3*WARP*WARP];
	
	float M, S, j;

	Sum_partials_nonregular( &M, &S, &j, d_input, s_input, size);
	
	if(d_pp[2]>0){
		Merge(&M, &S, &j, d_pp[0], d_pp[1], d_pp[2]);
	}
	
	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		d_MSD[0] = M / j;
		d_MSD[1] = sqrt(S / j);
		d_MSD[2] = j;
		d_pp[0] = M;
		d_pp[1] = S;
		d_pp[2] = j;
	}
}


__global__ void MSD_GPU_LA_ALL(float const* __restrict__ d_input, float *d_output, float *d_output_taps, int y_steps, int nTaps, int nTimesamples, int offset) {
	__shared__ float s_input[3*PD_NTHREADS];
	__shared__ float s_base[3*PD_NTHREADS];
	
	// MSD variables
	float M, S, j;
	float M_b, S_b, j_b;
	// FIR variables
	int d, gpos, spos, local_id;
	ushort EpT, limit;
	float2 ftemp1, ftemp2, ftemp3;
	float Bw[2];
	
	EpT = 2*PD_NTHREADS-nTaps+4;
	limit = blockDim.x - (nTaps>>2) - 1;

	// First y coordinate is separated
	//-------------------> FIR
	spos = blockIdx.x*(EpT) + 2*threadIdx.x;
	gpos = blockIdx.y*y_steps*nTimesamples + spos;
	Bw[0]=0; Bw[1]=0; j=0; j_b=0;
	if( (spos+4)<(nTimesamples-offset) ){
		// loading data for FIR filter. Each thread calculates two samples
		ftemp1.x= __ldg(&d_input[gpos]);	
		ftemp1.y= __ldg(&d_input[gpos+1]);
		ftemp2.x= __ldg(&d_input[gpos+2]);
		ftemp2.y= __ldg(&d_input[gpos+3]);
		ftemp3.x= __ldg(&d_input[gpos+4]);
		
		// Calculate FIR of 4 taps
		Bw[0]=ftemp1.x + ftemp1.y + ftemp2.x + ftemp2.y;
		Bw[1]=ftemp1.y + ftemp2.x + ftemp2.y + ftemp3.x;
		
		// Initialization of MSD variables for non-processed StrDev
		Initiate( &M_b, &S_b, &j_b, ftemp1.x );
		// First addition (second actually, but first done this way) non-processed StrDev
		Add_one( &M_b, &S_b, &j_b, ftemp1.y );
	}
	
	s_input[2*threadIdx.x] = Bw[0];
	s_input[2*threadIdx.x+1] = Bw[1];
	
	__syncthreads();
	
	// Calculating FIR up to nTaps
	for(d=4; d<nTaps; d=d+4){
		local_id = threadIdx.x+(d>>1);
		if( local_id<=limit ){
			Bw[0] = Bw[0] + s_input[2*local_id]; Bw[1] = Bw[1] + s_input[2*local_id+1];
		}
	}
	
	// Note: threads with local_id<0 which have wrong result create sums as well but are removed from final results later
	//       same is for base values as these would be included twice. First time here and next time in threadblock next to it
	//       this is due to halo needed for FIR filter	
	Initiate( &M, &S, &j, Bw[0] ); // Initialization of MSD variables for processed StrDev
	Add_one( &M, &S, &j, Bw[1] ); // First addition (second actually, but first done this way) processed StrDev
	
	
	// Rest of the iteration in y direction	
	for (int yf = 1; yf < y_steps; yf++) {
		__syncthreads();
		//-------------------> FIR
		spos = blockIdx.x*(EpT) + 2*threadIdx.x;
		gpos = blockIdx.y*y_steps*nTimesamples + yf*nTimesamples + spos;
		Bw[0]=0; Bw[1]=0;
		if( (spos+4)<(nTimesamples-offset) ){
			ftemp1.x= __ldg(&d_input[gpos]);	
			ftemp1.y= __ldg(&d_input[gpos+1]);
			ftemp2.x= __ldg(&d_input[gpos+2]);
			ftemp2.y= __ldg(&d_input[gpos+3]);
			ftemp3.x= __ldg(&d_input[gpos+4]);

			Bw[0]=ftemp1.x + ftemp1.y + ftemp2.x + ftemp2.y;
			Bw[1]=ftemp1.y + ftemp2.x + ftemp2.y + ftemp3.x;
			
			Add_one( &M_b, &S_b, &j_b, ftemp1.x );
			Add_one( &M_b, &S_b, &j_b, ftemp1.y );
		}
		
		s_input[2*threadIdx.x] = Bw[0];
		s_input[2*threadIdx.x+1] = Bw[1];
	
		__syncthreads();
	
		for(d=4; d<nTaps; d=d+4){	
			local_id = threadIdx.x+(d>>1);
			if( local_id<=limit ){
				Bw[0] = Bw[0] + s_input[2*local_id]; Bw[1] = Bw[1] + s_input[2*local_id+1];
			}
		}
		
		Add_one( &M, &S, &j, Bw[0] );
		Add_one( &M, &S, &j, Bw[1] );
	}
	
	__syncthreads();
	
	s_input[threadIdx.x] = 0;
	s_input[blockDim.x + threadIdx.x] = 0;
	s_input[2*blockDim.x + threadIdx.x] = 0;
	
	s_base[threadIdx.x] = 0;
	s_base[blockDim.x + threadIdx.x] = 0;
	s_base[2*blockDim.x + threadIdx.x] = 0;
	
	__syncthreads();

	spos=blockIdx.x*(EpT) + 2*threadIdx.x;	
	if( local_id<=limit ) {
		// Note: ommited number of samples in the last trailing threadblocks is due to -nTaps which is here. 
		//       Missing data should be contained in local_id. Thus this code is missing some time sample even it it does not need to. 
		//       When removed it produces different number of added time samples in j and j_b which is wierd
		if( spos<(nTimesamples-offset-nTaps) ) { // -nTaps might not be necessary
			s_input[local_id] = M;
			s_input[blockDim.x + local_id] = S;
			s_input[2*blockDim.x + local_id] = j;
			
			s_base[local_id] = M_b;
			s_base[blockDim.x + local_id] = S_b;
			s_base[2*blockDim.x + local_id] = j_b;
		}

	}
	__syncthreads();
	
	//------------------------------------------------------------------------------------
	//---------> StrDev of processed input
	Reduce_SM( &M, &S, &j, s_input );
	Reduce_WARP( &M, &S, &j);
	
	//------------------------------------------------------------------------------------
	//---------> StrDev of unprocessed input
	Reduce_SM( &M_b, &S_b, &j_b, s_base );
	Reduce_WARP( &M_b, &S_b, &j_b);

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


__global__ void MSD_GPU_LA_ALL_Nth(float const* __restrict__ d_input, float const* __restrict   d_bv_in, float *d_output, float *d_output_taps, int y_steps, int nTaps, int nTimesamples, int offset) {
	__shared__ float s_input[3*PD_NTHREADS];
	__shared__ float s_base[3*PD_NTHREADS];
	
	// MSD variables
	float M, S, j;
	float M_b, S_b, j_b;
	// FIR variables
	int d, gpos, spos, local_id;
	ushort EpT, limit;
	float2 ftemp1, ftemp2, ftemp3;
	float Bw[2];
	
	EpT = 2*PD_NTHREADS-nTaps+4;
	limit = blockDim.x - (nTaps>>2) - 1;

	// First y coordinate is separated
	//-------------------> FIR
	spos = blockIdx.x*(EpT) + 2*threadIdx.x;
	gpos = blockIdx.y*y_steps*nTimesamples + spos;
	Bw[0]=0; Bw[1]=0; j=0; j_b=0;
	if( (spos+4)<(nTimesamples-offset) ){
		// loading data for FIR filter. Each thread calculates two samples
		ftemp1.x= __ldg(&d_input[gpos]);	
		ftemp1.y= __ldg(&d_input[gpos+1]);
		ftemp2.x= __ldg(&d_input[gpos+2]);
		ftemp2.y= __ldg(&d_input[gpos+3]);
		ftemp3.x= __ldg(&d_input[gpos+4]);
		
		// Calculate FIR of 4 taps
		Bw[0]=ftemp1.x + ftemp1.y + ftemp2.x + ftemp2.y;
		Bw[1]=ftemp1.y + ftemp2.x + ftemp2.y + ftemp3.x;
		
		Initiate( &M_b, &S_b, &j_b, __ldg(&d_bv_in[gpos]) );
		Add_one( &M_b, &S_b, &j_b, __ldg(&d_bv_in[gpos+1]) );
	}
	
	s_input[2*threadIdx.x] = Bw[0];
	s_input[2*threadIdx.x+1] = Bw[1];
	
	__syncthreads();
	
	// Calculating FIR up to nTaps
	for(d=4; d<nTaps; d=d+4){
		local_id = threadIdx.x+(d>>1);
		if( local_id<=limit ){
			Bw[0] = Bw[0] + s_input[2*local_id]; Bw[1] = Bw[1] + s_input[2*local_id+1];
		}
	}
	
	// Note: threads with local_id<0 which have wrong result create sums as well but are removed from final results later
	//       same is for base values as these would be included twice. First time here and next time in threadblock next to it
	//       this is due to halo needed for FIR filter
	Initiate( &M, &S, &j, __ldg(&d_bv_in[gpos]) + Bw[0] );
	Add_one( &M, &S, &j, __ldg(&d_bv_in[gpos+1]) + Bw[1] );
	
	
	// Rest of the iteration in y direction	
	for (int yf = 1; yf < y_steps; yf++) {
		__syncthreads();
		//-------------------> FIR
		spos = blockIdx.x*(EpT) + 2*threadIdx.x;
		gpos = blockIdx.y*y_steps*nTimesamples + yf*nTimesamples + spos;
		Bw[0]=0; Bw[1]=0;
		if( (spos+4)<(nTimesamples-offset) ){
			ftemp1.x= __ldg(&d_input[gpos]);	
			ftemp1.y= __ldg(&d_input[gpos+1]);
			ftemp2.x= __ldg(&d_input[gpos+2]);
			ftemp2.y= __ldg(&d_input[gpos+3]);
			ftemp3.x= __ldg(&d_input[gpos+4]);

			Bw[0]=ftemp1.x + ftemp1.y + ftemp2.x + ftemp2.y;
			Bw[1]=ftemp1.y + ftemp2.x + ftemp2.y + ftemp3.x;
			
			Add_one( &M_b, &S_b, &j_b, __ldg(&d_bv_in[gpos]) );
			Add_one( &M_b, &S_b, &j_b, __ldg(&d_bv_in[gpos+1]) );
		}
		
		s_input[2*threadIdx.x] = Bw[0];
		s_input[2*threadIdx.x+1] = Bw[1];
	
		__syncthreads();
	
		for(d=4; d<nTaps; d=d+4){	
			local_id = threadIdx.x+(d>>1);
			if( local_id<=limit ){
				Bw[0] = Bw[0] + s_input[2*local_id]; Bw[1] = Bw[1] + s_input[2*local_id+1];
			}
		}
		
		Add_one( &M, &S, &j, __ldg(&d_bv_in[gpos]) + Bw[0] );
		Add_one( &M, &S, &j, __ldg(&d_bv_in[gpos+1]) + Bw[1] );
	}
	
	__syncthreads();
	
	s_input[threadIdx.x] = 0;
	s_input[blockDim.x + threadIdx.x] = 0;
	s_input[2*blockDim.x + threadIdx.x] = 0;
	
	s_base[threadIdx.x] = 0;
	s_base[blockDim.x + threadIdx.x] = 0;
	s_base[2*blockDim.x + threadIdx.x] = 0;
	
	__syncthreads();

	spos=blockIdx.x*(EpT) + 2*threadIdx.x;	
	if( local_id<=limit ) {		
		// Note: ommited number of samples in the last trailing threadblocks is due to -nTaps which is here. 
		//       Missing data should be contained in local_id. Thus this code is missing some time sample even it it does not need to. 
		//       When removed it produces different number of added time samples in j and j_b which is wierd
		if( spos<(nTimesamples-offset-nTaps) ) { // -nTaps might not be necessary
			s_input[local_id] = M;
			s_input[blockDim.x + local_id] = S;
			s_input[2*blockDim.x + local_id] = j;
			
			s_base[local_id] = M_b;
			s_base[blockDim.x + local_id] = S_b;
			s_base[2*blockDim.x + local_id] = j_b;
		}

	}
	__syncthreads();
	
	//------------------------------------------------------------------------------------
	//---------> StrDev of processed input
	Reduce_SM( &M, &S, &j, s_input );
	Reduce_WARP( &M, &S, &j);
	
	//------------------------------------------------------------------------------------
	//---------> StrDev of unprocessed input
	Reduce_SM( &M_b, &S_b, &j_b, s_base );
	Reduce_WARP( &M_b, &S_b, &j_b);

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


__global__ void MSD_GPU_final_create_LA(float *d_input, float *d_output, float *d_MSD_base, int nTaps, int size) {
	__shared__ float s_input[3*WARP*WARP];

	float M, S, j;

	Sum_partials_regular( &M, &S, &j, d_input, s_input, size);
	
	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		d_output[0] = d_MSD_base[0];
		d_output[1] = d_MSD_base[1];
		d_output[2] = (sqrt(S / j) - d_MSD_base[1])/( (float) (nTaps-1));
	}
}


__global__ void MSD_GPU_final_create_LA_Nth(float *d_input, float *d_output, float *d_MSD_base, float *d_MSD_DIT, int nTaps, int size, int DIT_value) {
	__shared__ float s_input[3*WARP*WARP];

	float M, S, j;

	Sum_partials_regular( &M, &S, &j, d_input, s_input, size);

	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0) {
		d_output[0] = d_MSD_base[0];
		d_output[1] = d_MSD_base[1];
		d_output[2] = (sqrt(S / j) - d_MSD_base[1])/( (float) nTaps);
		d_output[3] = d_MSD_DIT[0]*DIT_value; 
	}
}


#endif
