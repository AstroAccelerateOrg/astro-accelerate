#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_log.hpp"
#include "aa_params.hpp"
#include "aa_device_cuda_deprecated_wrappers.cuh"

namespace astroaccelerate {

#define MEAN	127.5f
#define CUT 	2.0f
#define R_CUT	4.0f
#define ITER 	20
#define ACC	0.001f

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


__device__ __inline__ void MSD_block(unsigned short *d_input, int *x_steps, int *nSamples, float *s_par_MSD, int *s_par_nElements, float *mean, float *stdev){	
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


__device__ __inline__ void MSD_block_outlier_rejection(unsigned short *d_input, int *x_steps, int *nSamples, float *sigma, float *s_par_MSD, int *s_par_nElements, float *mean, float *stdev){	
	float M, S, ftemp;
	int j;
	
	float limit_down = (*mean) - (*sigma)*(*stdev);
	float limit_up = (*mean) + (*sigma)*(*stdev);

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


__global__ void zero_dm_outliers_kernel_channels(unsigned short *d_input, int nchans, float outlier_sigma, int nbits, float *normalization_factor) {
	__shared__ float s_par_MSD[256];
	__shared__ int s_par_nElements[128];
	float mean, stdev, old_mean;
	unsigned long int gpos = blockIdx.x*nchans;
	int x_steps = (nchans+blockDim.x-1)/blockDim.x;
	mean =0; stdev=0; old_mean=0;

	MSD_block(&d_input[gpos], &x_steps, &nchans, s_par_MSD, s_par_nElements, &mean, &stdev);
	
	for(int f = 0; f < 5; f++){
		old_mean = mean;
		MSD_block_outlier_rejection(&d_input[gpos], &x_steps, &nchans, &outlier_sigma, s_par_MSD, s_par_nElements, &mean, &stdev);
		if(old_mean-mean < ACC) break;
	}
	
	for(int c = 0; c < x_steps; c++){
		if ((c*blockDim.x + threadIdx.x) < nchans) {
			unsigned short value;
			float result = (float)d_input[blockIdx.x*nchans + c*blockDim.x + threadIdx.x] - mean + normalization_factor[c*blockDim.x + threadIdx.x];
			if( nbits == 4 ){
				if(result<0) value = 0;
				else if(result>15) value = 15;
				else value = (unsigned short) result;
			}
			else if( nbits == 8 ){
				if(result<0) value = 0;
				else if(result>255) value = 255;
				else value = (unsigned short) result;
			}
			else if( nbits == 16 ){
				if(result<0) value = 0;
				else if(result>65535) value = 65535;
				else value = (unsigned short) result;
			}
			d_input[blockIdx.x*nchans + c*blockDim.x + threadIdx.x] = value;
		}
	}
}

/**
* \brief Dristribution based RFI removal.
* \todo Needs cleaning and optimizing (WA 07/08/18)
*/
__global__ void zero_dm_outliers_kernel_one(unsigned short *d_input, int nchans, int nsamp) {
	
	int t  = blockIdx.x * blockDim.x + threadIdx.x;
	
	int count = 0;
	int iters = 0;
	
	float stdev = 1000000.0f;
	float mean = MEAN;
	float mean_last = 0.0f;
	float sum = 0.0f;
	float sum_squares = 0.0f;
	float cutoff = (CUT * stdev);
	
	__shared__ float g_mean[1024];
	__shared__ float g_stdev[1024];
	
	// Calculation of MSD until convergence
	while(abs(mean - mean_last) > ACC) {
		sum = 0.0f;
		sum_squares = 0.0f;
		count = 0;
		
		for(int c = 0; c < nchans; c++) {
			float data=(float)d_input[t*nchans + c];
			if(data < (mean + cutoff) && data > (mean - cutoff) ) {
				sum+=data;
				sum_squares+=(data*data);
				count++;
			}
		}
		mean_last = mean;
		mean = (sum/(float)count);
		sum_squares = ((sum_squares / count) - (mean * mean));
		stdev = sqrt(sum_squares);
		cutoff = (CUT * stdev);
		
		iters++;
		if(iters > ITER) break;
	}
	//-------------------------------------------------------------<
	
	// Set channel values to default mean or calculated mean
	if(count == 0 || iters > ITER || mean == 0.0f || stdev == 0.0f) {
		for(int c = 0; c < nchans; c++) {
			d_input[t*nchans + c] = MEAN;
		}
		g_mean[threadIdx.x] = mean = MEAN;
		g_stdev[threadIdx.x] = stdev = 0.0f;
	} 
	else {
		g_mean[threadIdx.x] = mean;
		g_stdev[threadIdx.x] = stdev;
	}
	
	__syncthreads();
	
	float 	mean_of_mean = 0.0f;
	float 	stdev_of_mean = 0.0f;
	float 	m_cutoff = 0.0f;
	
	sum_squares = 0.0f;
	
	for(int i = 0; i<blockDim.x; i++) {
		mean_of_mean += g_mean[i];
		sum_squares   += (g_mean[i]* g_mean[i]);
	}
	
	mean_of_mean /= blockDim.x;
	sum_squares = ((sum_squares / blockDim.x) - (mean_of_mean * mean_of_mean));
	
	stdev_of_mean = sqrt(sum_squares);
	
	m_cutoff = (3.0*stdev_of_mean);
	
	float 	mean_of_stdev = 0.0f;
	float 	stdev_of_stdev = 0.0f;
	float 	s_cutoff = 0.0f;
	
	sum_squares = 0.0f;
	
	for(int i = 0; i<blockDim.x; i++) {
		mean_of_stdev += g_stdev[i];
		sum_squares   += (g_stdev[i]* g_stdev[i]);
	}
	
	mean_of_stdev /= blockDim.x;
	sum_squares = ((sum_squares / blockDim.x) - (mean_of_stdev * mean_of_stdev));
	
	stdev_of_stdev = sqrt(sum_squares);
	
	s_cutoff = (3.0*stdev_of_stdev);
	
	
	if((g_mean[threadIdx.x] - mean_of_mean) > m_cutoff ||  (g_mean[threadIdx.x] - mean_of_mean) < -m_cutoff) {
		for(int c = 0; c < nchans; c++) {
			d_input[t*nchans + c] = MEAN;
		}
	} 
	else if((g_stdev[threadIdx.x] - mean_of_stdev) > s_cutoff ||  (g_stdev[threadIdx.x] - mean_of_stdev) < -s_cutoff) {
		for(int c = 0; c < nchans; c++) {
			d_input[t*nchans + c] = MEAN;
		}
	} 
	else {
		for(int c = 0; c < nchans; c++) {
			if((d_input[t*nchans + c]-mean < R_CUT*stdev) && (d_input[t*nchans + c]-mean > - R_CUT*stdev)) {
				d_input[t*nchans + c] = (unsigned short)((float)d_input[t*nchans + c]-(float)mean+MEAN);
			} 
			else {
				d_input[t*nchans + c] = MEAN;
			}
		}
	}
}


__global__ void zero_dm_outliers_kernel_two(unsigned short *d_input, int nchans, int nsamp) {
	int count = 0;
	int iters = 0;

	float stdev = 1000000.0f;
	float mean = MEAN;
	float mean_last = 0.0f;
	float sum = 0.0f;
	float sum_squares = 0.0f;
	float cutoff = (CUT * stdev);


	__shared__ float g_mean[1024];
	__shared__ float g_stdev[1024];


	int c = blockIdx.x * blockDim.x + threadIdx.x;

	count = 0;
	iters = 0;

	stdev = 1000000.0f;
	mean = MEAN;
	mean_last = 0.0f;
	cutoff = (CUT * stdev);

	while(abs(mean - mean_last) > ACC) {
		sum = 0.0f;
		sum_squares = 0.0f;
		count = 0;

		for(int t = 0; t < nsamp; t++) {
			float data=(float)d_input[t*nchans + c];
			if(data < (mean + cutoff) && data > (mean - cutoff) ) {
				sum+=data;
				sum_squares+=(data*data);
				count++;
			}
		}
		mean_last = mean;
		mean = (sum/(float)count);
		sum_squares = ((sum_squares / count) - (mean * mean));
		stdev = sqrt(sum_squares);
		cutoff = (CUT * stdev);

		iters++;
		if(iters > ITER) break;
	}

	if(count == 0 || iters > ITER || mean == 0.0f || stdev == 0.0f) {
		for(int t = 0; t < nsamp; t++) {
			d_input[t*nchans + c] = MEAN;
		}
		g_mean[threadIdx.x] = mean = MEAN;
		g_stdev[threadIdx.x] = stdev = 0.0f;
	} 
	else {
		g_mean[threadIdx.x] = mean;
		g_stdev[threadIdx.x] = stdev;
	}

	__syncthreads();

	float 	mean_of_mean = 0.0f;
	float 	stdev_of_mean = 0.0f;
	float 	m_cutoff = 0.0f;

	sum_squares = 0.0f;

	for(int i = 0; i<blockDim.x; i++) {
		mean_of_mean += g_mean[i];
		sum_squares   += (g_mean[i]* g_mean[i]);
	}

	mean_of_mean /= blockDim.x;
	sum_squares = ((sum_squares / blockDim.x) - (mean_of_mean * mean_of_mean));

	stdev_of_mean = sqrt(sum_squares);

	m_cutoff = (3.0*stdev_of_mean);

	float 	mean_of_stdev = 0.0f;
	float 	stdev_of_stdev = 0.0f;
	float 	s_cutoff = 0.0f;

	sum_squares = 0.0f;

	for(int i = 0; i<blockDim.x; i++) {
		mean_of_stdev += g_stdev[i];
		sum_squares   += (g_stdev[i]* g_stdev[i]);
	}

	mean_of_stdev /= blockDim.x;
	sum_squares = ((sum_squares / blockDim.x) - (mean_of_stdev * mean_of_stdev));

	stdev_of_stdev = sqrt(sum_squares);

	s_cutoff = (3.0*stdev_of_stdev);

	if((g_mean[threadIdx.x] - mean_of_mean) > m_cutoff ||  (g_mean[threadIdx.x] - mean_of_mean) < -m_cutoff) {
		for(int t = 0; t < nsamp; t++) {
			d_input[t*nchans + c] = MEAN;
		}
	} 
	else if((g_stdev[threadIdx.x] - mean_of_stdev) > s_cutoff ||  (g_stdev[threadIdx.x] - mean_of_stdev) < -s_cutoff) {
		for(int t = 0; t < nsamp; t++) {
			d_input[t*nchans + c] = MEAN;
		}
	}
	else {
		for(int t = 0; t < nsamp; t++) {
			if((d_input[t*nchans + c]-mean < R_CUT*stdev) && (d_input[t*nchans + c]-mean > - R_CUT*stdev)) {
				d_input[t*nchans + c] = (unsigned short)((float)d_input[t*nchans + c]-(float)mean+MEAN);
			} 
			else {
				d_input[t*nchans + c] = MEAN;
			}
		}
	}
}


//---------------------------- WRAPPERS --------------------------------

void call_kernel_zero_dm_outliers_kernel_channels(
	const dim3 &grid_size, 
	const dim3 &block_size, 
	const int &smem, 
	const cudaStream_t &stream, 
	unsigned short *const d_input, 
	const int &nchans, 
	const float &outlier_sigma, 
	const int &nbits, 
	float *normalization_factor) {
		if (block_size.x<=0 || block_size.y<=0 || block_size.z<=0 || grid_size.x<=0 || grid_size.y<=0 || grid_size.z<=0){
			LOG(log_level::error, "Zero DM kernel is configured incorrectly. grid_size=[" + std::to_string(grid_size.x) + "; " + std::to_string(grid_size.y) + "; " + std::to_string(grid_size.z) + "] block_size=[" + std::to_string(block_size.x) + "; " + std::to_string(block_size.y) + "; " + std::to_string(block_size.z) + "]");
		}
		else {
			zero_dm_outliers_kernel_channels<<< grid_size, block_size, smem, stream >>>(d_input, nchans, outlier_sigma, nbits, normalization_factor);
			cudaDeviceSynchronize();
		}
}

/** \brief Kernel wrapper function for zero_dm_outliers_kernel_one kernel function. */
void call_kernel_zero_dm_outliers_kernel_one(const dim3 &block_size, const dim3 &grid_size,
unsigned short *const d_input, const int &nchans, const int &nsamp) {
zero_dm_outliers_kernel_one<<<block_size, grid_size>>>(d_input, nchans, nsamp);
}

/** \brief Kernel wrapper function for zero_dm_outliers_kernel_two kernel function. */
void call_kernel_zero_dm_outliers_kernel_two(const dim3 &block_size, const dim3 &grid_size,
unsigned short *const d_input, const int &nchans, const int &nsamp) {
zero_dm_outliers_kernel_two<<<block_size, grid_size>>>(d_input, nchans, nsamp);
}

} //namespace astroaccelerate
