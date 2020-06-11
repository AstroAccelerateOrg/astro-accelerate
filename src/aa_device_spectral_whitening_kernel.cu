#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "aa_params.hpp"
#include "aa_device_cuda_deprecated_wrappers.cuh"

namespace astroaccelerate {

__device__ __inline__ void Initiate(float *M, float *S, int *j, float element){
	*M = element;
	*S = 0;
	*j = 1;
}

__device__ __inline__ void Add_one(float *M, float *S, int *j, float element){
	float ftemp;
	*j = (*j) + 1;
	double r_j = (double) (*j);
	*M = (*M) + element;
	ftemp = ( r_j*element - (*M) );
	*S = (*S) + 1.0f / ( r_j*( r_j - 1.0f ) )*ftemp*ftemp;
}

__device__ __inline__ void Merge(float *A_M, float *A_S, int *A_j, float B_M, float B_S, int B_j){
	float ftemp;
	double r_B_j = (double) B_j;
	double r_A_j = (double) (*A_j);
	
	ftemp = ( r_B_j / r_A_j)*(*A_M) - B_M;
	(*A_S) = (*A_S) + B_S + ( r_A_j/( r_B_j*(r_A_j + r_B_j) ) )*ftemp*ftemp;
	(*A_M) = (*A_M) + B_M;
	(*A_j) = (*A_j) + B_j;
}

__device__ __inline__ void Reduce_SM(float *M, float *S, int *j, float *s_par_MSD, int *s_par_nElements){
	int jv;
	
	(*M)=s_par_MSD[threadIdx.x];
	(*S)=s_par_MSD[blockDim.x + threadIdx.x];
	(*j)=s_par_nElements[threadIdx.x];
	
	for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1) {
		if (threadIdx.x < i) {
			jv = s_par_nElements[i + threadIdx.x];
			if( jv!=0){
				if( (*j)==0 ){
					(*M) = s_par_MSD[i + threadIdx.x];
					(*S) = s_par_MSD[blockDim.x + i + threadIdx.x];
					(*j) = jv;
				}
				else {
					Merge(M, S, j, s_par_MSD[i + threadIdx.x], s_par_MSD[blockDim.x + i + threadIdx.x], jv);
				}
			}
			
			s_par_MSD[threadIdx.x] = (*M);
			s_par_MSD[blockDim.x + threadIdx.x] = (*S);
			s_par_nElements[threadIdx.x] = (*j);
		}
		__syncthreads();
	}
}

__device__ __inline__ void Reduce_WARP(float *M, float *S, int *j){
	float B_M, B_S;
	int B_j;
	
	for (int q = HALF_WARP; q > 0; q = q >> 1) {
		B_M = aa_shfl_down(AA_ASSUME_MASK, (*M), q);
		B_S = aa_shfl_down(AA_ASSUME_MASK, (*S), q);
		B_j = aa_shfl_down(AA_ASSUME_MASK, (*j), q);
		
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


__device__ __inline__ void MSD_block(float *d_input, int *x_steps, unsigned long int *nSamples, float *s_par_MSD, int *s_par_nElements, float *mean, float *stdev){	
	float M, S, ftemp;
	int j;	

	M=0; S=0; j=0;
	int spos = threadIdx.x;
	if( spos < (*nSamples) ){
		ftemp = (float) d_input[spos];
		Initiate( &M, &S, &j, ftemp);
		
		spos = spos + blockDim.x;
		for (int xf = 1; xf < (*x_steps); xf++) {
			if( spos < (*nSamples) ){
				ftemp = (float) d_input[spos];
				Add_one( &M, &S, &j, ftemp);
				spos = spos + blockDim.x;
			}
		}
	}
	
	s_par_MSD[threadIdx.x] = M;
	s_par_MSD[128 + threadIdx.x] = S;
	s_par_nElements[threadIdx.x] = j;
	
	__syncthreads();
	
	Reduce_SM( &M, &S, &j, s_par_MSD, s_par_nElements );
	Reduce_WARP( &M, &S, &j);
	
	if(threadIdx.x == 0){
		(*mean) = M / (double) j;
		(*stdev) = sqrt(S / (double) j);
		s_par_MSD[0] = (*mean);
		s_par_MSD[1] = (*stdev);
	}
	
	__syncthreads();
	
	(*mean) = s_par_MSD[0];
	(*stdev) = s_par_MSD[1];
	
	__syncthreads();
}


__device__ __inline__ void MSD_block_outlier_rejection(float *d_input, int *x_steps, unsigned long int *nSamples, float *s_par_MSD, int *s_par_nElements, float *mean, float *stdev){	
	float M, S, ftemp;
	int j;
	
	float sigma = 3.0;
	float limit_down = (*mean) - sigma*(*stdev);
	float limit_up = (*mean) + sigma*(*stdev);

	M=0;	S=0;	j=0;
	int spos = threadIdx.x;
	if( spos < (*nSamples) ){
		for (int xf = 0; xf < (*x_steps); xf++) {
			if( spos < (*nSamples) ){
				ftemp = (float) d_input[spos];
				if( (ftemp>limit_down) && (ftemp < limit_up) ){
					if(j==0){
						Initiate( &M, &S, &j, ftemp);
					}
					else{
						Add_one( &M, &S, &j, ftemp);
					}			
				}
				spos = spos + blockDim.x;
			}
		}	
	}

	s_par_MSD[threadIdx.x] = M;
	s_par_MSD[blockDim.x + threadIdx.x] = S;
	s_par_nElements[threadIdx.x] = j;
	
	__syncthreads();
	
	Reduce_SM( &M, &S, &j, s_par_MSD, s_par_nElements );
	Reduce_WARP( &M, &S, &j);
	
	if(threadIdx.x == 0){
		(*mean) = M / (double) j;
		(*stdev) = sqrt(S / (double) j);
		s_par_MSD[0] = (*mean);
		s_par_MSD[1] = (*stdev);
	}
	
	__syncthreads();
	
	(*mean) = s_par_MSD[0];
	(*stdev) = s_par_MSD[1];
	
	__syncthreads();
}


__global__ void GPU_kernel_spectrum_whitening_SGP1(float *d_input, unsigned long int nSamples) {
	__shared__ float s_par_MSD[256];
	__shared__ int s_par_nElements[128];

	float mean, stdev;
	unsigned long int gpos = blockIdx.y*nSamples;
	unsigned long int spos = blockIdx.x*blockDim.x;
	int x_steps = 1;
	
	MSD_block(&d_input[gpos + spos], &x_steps, &nSamples, s_par_MSD, s_par_nElements, &mean, &stdev);
	
	MSD_block_outlier_rejection(&d_input[gpos + spos], &x_steps, &nSamples, s_par_MSD, s_par_nElements, &mean, &stdev);
	
	d_input[gpos + spos + threadIdx.x] = (d_input[gpos + spos + threadIdx.x] - mean)/stdev;
}

//---------------------------- WRAPPERS --------------------------

void call_kernel_spectrum_whitening_SGP1(
	const dim3 &block_size, 
	const dim3 &grid_size, 
	const int &smem_bytes, 
	const cudaStream_t &stream, 
	float *d_input, 
	unsigned long int nSamples
) {	
	GPU_kernel_spectrum_whitening_SGP1<<< grid_size, block_size, smem_bytes, stream >>>(d_input, nSamples);
}

} //namespace astroaccelerate

