#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include "aa_params.hpp"
#include "aa_device_cuda_deprecated_wrappers.cuh"

namespace astroaccelerate {

class Sort_Params {
public:
	static const int warp = 32;
	// search_exp = log_2(N) + 1
};

class Sort_64 : public Sort_Params {
public:
	static const int nThreads = 32;
	static const int search_exp = 7;
};

class Sort_128 : public Sort_Params {
public:
	static const int nThreads = 64;
	static const int search_exp = 8;
};

class Sort_256 : public Sort_Params {
public:
	static const int nThreads = 128;
	static const int search_exp = 9;
};

class Sort_512 : public Sort_Params {
public:
	static const int nThreads = 256;
	static const int search_exp = 10;
};

class Sort_1024 : public Sort_Params {
public:
	static const int nThreads = 512;
	static const int search_exp = 11;
};

class Sort_2048 : public Sort_Params {
public:
	static const int nThreads = 1024;
	static const int search_exp = 12;
};


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

//--------------- Bitonic sort device functions -----------------

template<typename intype>
__device__ __inline__ void butterfly_exchange(intype *s_values, int butterfly_size, int sort_direction){
	int block_size      = (1<<butterfly_size);
	int block_size_half = (1<<(butterfly_size-1));
	int local_block     = ((threadIdx.x)>>(butterfly_size-1));
	int local_thread    = ((threadIdx.x)&(block_size_half-1));
	
	int elA, elB;
	elA = local_block*block_size + local_thread;
	elB = local_block*block_size + local_thread + (block_size>>1);
	
	if(sort_direction==0 && s_values[elA] > s_values[elB]){ // ascending
		intype swap;
		swap = s_values[elA];
		s_values[elA] = s_values[elB];
		s_values[elB] = swap;
	}
	if(sort_direction==1 && s_values[elA] < s_values[elB]){ // descending
		intype swap;
		swap = s_values[elA];
		s_values[elA] = s_values[elB];
		s_values[elB] = swap;
	}
}

template<class const_params, typename intype>
__device__ __inline__ void bitonic_sort_block(intype *s_values){
	// major cycle
	//     Consists of series of incrementally larger butterfly operations 
	for(int majorc = 1; majorc<const_params::search_exp; majorc++){
		int sort_dir_flag = (threadIdx.x>>(majorc-1))&1;
		// minor cycle
		//     Consists of series of butterfly operations of smaller and smaller size
		//     until size of two
		for(int minorc = majorc; minorc>0; minorc--){
			butterfly_exchange(s_values, minorc, sort_dir_flag);
			__syncthreads();
		}
	}
}

//--------------- Bitonic sort device functions -----------------



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


__global__ void GPU_kernel_segmented_MSD(float *d_segmented_MSD, float2 *d_input, int *d_segment_sizes, size_t nSamples, int nSegments){
	float M, S, ftemp;
	int j, spos;
	size_t gpos;
	
	int start_pos = d_segment_sizes[blockIdx.x];
	int end_pos   = d_segment_sizes[blockIdx.x + 1];
	int active_range = end_pos - start_pos;
	
	
	M=0; S=0; j=0;
	spos = threadIdx.x;
	
	//first iteration
	gpos = blockIdx.y*nSamples + start_pos + spos;
	if(spos < active_range){
		ftemp = d_input[gpos].x*d_input[gpos].x + d_input[gpos].y*d_input[gpos].y;
		if(blockIdx.x==0 && threadIdx.x==0){
			ftemp = 0;
		}
		Initiate( &M, &S, &j, ftemp);
	}
	spos = spos + 32;
	
	__syncthreads();
	
	for(int i = 0; i < 7; i++){
		if(spos < active_range){
			gpos = blockIdx.y*nSamples + start_pos + spos;
			ftemp = d_input[gpos].x*d_input[gpos].x + d_input[gpos].y*d_input[gpos].y;
			Add_one( &M, &S, &j, ftemp);
		}
		spos = spos + 32;
	}
	
	__syncthreads();
	
	Reduce_WARP( &M, &S, &j);
	
	if(threadIdx.x == 0){
		d_segmented_MSD[blockIdx.y*nSegments  + blockIdx.x] = M / (double) j;
		//d_segmented_MSD[blockIdx.y*nSegments  + 2*blockIdx.x + 1] = sqrt(S / (double) j);
	}
}

template<class const_params>
__global__ void GPU_kernel_segmented_median(float *d_segmented_MSD, float2 *d_input, int *d_segment_sizes, size_t nSamples, int nSegments){
	__shared__ float s_values[2*const_params::nThreads];
	s_values[threadIdx.x] = 0;
	s_values[threadIdx.x + const_params::nThreads] = 0;
	
	int start_pos = d_segment_sizes[blockIdx.x];
	int end_pos   = d_segment_sizes[blockIdx.x + 1];
	int active_range = end_pos - start_pos;
	if(active_range <= 0 && active_range > const_params::nThreads*2) return;
	
	int local_pos;
	local_pos = start_pos + threadIdx.x;
	if( threadIdx.x < active_range && local_pos < nSamples) {
		float2 value = d_input[blockIdx.y*nSamples + local_pos];
		s_values[threadIdx.x] = value.x*value.x + value.y*value.y;
	}
	if( (threadIdx.x + const_params::nThreads) < active_range && (local_pos +  + const_params::nThreads) < nSamples ) {
		float2 value = d_input[blockIdx.y*nSamples + const_params::nThreads + local_pos];
		s_values[const_params::nThreads + threadIdx.x] = value.x*value.x + value.y*value.y;
	}
	__syncthreads();
	
	bitonic_sort_block<const_params>(s_values);
	
	int nZeros = const_params::nThreads*2 - active_range;
	int median_pos = nZeros + (active_range>>1);
	
	if( threadIdx.x == 0 ){
		if( (active_range&1)==0 ){ // even
			float median = (s_values[median_pos] + s_values[median_pos + 1])*0.5f;
			d_segmented_MSD[blockIdx.y*nSegments  + blockIdx.x] = median;
		}
		else {
			float median = s_values[median_pos];
			d_segmented_MSD[blockIdx.y*nSegments  + blockIdx.x] = median;
		}
	}
}


__global__ void GPU_kernel_spectrum_whitening_SGP2(float *d_segmented_MSD, float2 *d_input, int *d_segment_sizes, size_t nSamples, int nSegments){
	if(blockIdx.x==0){
		// First segment cannot calculate slope thus it is normalized by mean only
		int f0_pos = d_segment_sizes[blockIdx.x + 0];
		int fp1_pos = d_segment_sizes[blockIdx.x + 1];
		int current_range  = fp1_pos - f0_pos;
		float current_mean  = d_segmented_MSD[blockIdx.y*nSegments + blockIdx.x]*1.44269504088896;
		float norm = 1.0/sqrt(current_mean);
		size_t global_pos = blockIdx.y*nSamples + threadIdx.x;
		if(threadIdx.x<(current_range>>1)){
			d_input[global_pos].x *= norm;
			d_input[global_pos].y *= norm;
		}
		if(threadIdx.x==0){
			d_input[global_pos].x = 1.0;
			d_input[global_pos].y = 0.0;
		}
	}
	else {
		// This is presto deredning scheme
		int fm1_pos = d_segment_sizes[blockIdx.x - 1];
		int f0_pos  = d_segment_sizes[blockIdx.x + 0];
		int fp1_pos = d_segment_sizes[blockIdx.x + 1];
		int previous_range = f0_pos - fm1_pos;
		int current_range  = fp1_pos - f0_pos;
		float previous_mean = d_segmented_MSD[blockIdx.y*nSegments + blockIdx.x - 1]*1.44269504088896;
		float current_mean  = d_segmented_MSD[blockIdx.y*nSegments + blockIdx.x]*1.44269504088896;
		
		int range = ((previous_range + current_range)>>1);
		int i = range - threadIdx.x;
		float slope = (current_mean - previous_mean) / ((float) range);
		float norm  = 1.0/sqrt(previous_mean + slope*i);
		int local_pos = fm1_pos + (previous_range>>1) + i;
		
		if( i >= 0 && i<=range && local_pos < nSamples){
			size_t global_pos = blockIdx.y*nSamples + local_pos;
			d_input[global_pos].x *= norm;
			d_input[global_pos].y *= norm;
		}
		
		if(blockIdx.x == (gridDim.x - 2)){
			int local_pos = fm1_pos + (previous_range>>1) + range + threadIdx.x;
			if(local_pos < nSamples){
				size_t global_pos = blockIdx.y*nSamples + local_pos;
				d_input[global_pos].x *= norm;
				d_input[global_pos].y *= norm;
			}
		}
	}
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

void call_kernel_segmented_MSD(
	const  dim3 &block_size, 
	const  dim3 &grid_size, 
	const  int &smem_bytes, 
	const  cudaStream_t &stream, 
	float  *d_segmented_MSD, 
	float2 *d_input, 
	int    *d_segment_sizes,
	size_t nSamples,
	int    nSegments
) {	
	GPU_kernel_segmented_MSD<<< grid_size, block_size, smem_bytes, stream >>>(d_segmented_MSD, d_input, d_segment_sizes, nSamples, nSegments);
}

void call_kernel_segmented_median(
	const  dim3 &block_size, 
	const  dim3 &grid_size, 
	const  int &smem_bytes, 
	const  cudaStream_t &stream, 
	float  *d_segmented_MSD, 
	float2 *d_input, 
	int    *d_segment_sizes,
	size_t nSamples,
	int    nSegments
) {	
	GPU_kernel_segmented_median<Sort_256><<< grid_size, block_size, smem_bytes, stream >>>(d_segmented_MSD, d_input, d_segment_sizes, nSamples, nSegments);
}

void call_kernel_spectrum_whitening_SGP2(
	const  dim3 &block_size, 
	const  dim3 &grid_size, 
	const  int &smem_bytes, 
	const  cudaStream_t &stream, 
	float  *d_segmented_MSD, 
	float2 *d_input,
	int    *d_segment_sizes,
	size_t nSamples,
	int    nSegments
) {	
	GPU_kernel_spectrum_whitening_SGP2<<< grid_size, block_size, smem_bytes, stream >>>(d_segmented_MSD, d_input, d_segment_sizes, nSamples, nSegments);
}

} //namespace astroaccelerate





