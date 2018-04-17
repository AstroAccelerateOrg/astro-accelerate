#ifndef MSD_SHARED_KERNEL_FUNCTIONS_H_
#define MSD_SHARED_KERNEL_FUNCTIONS_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"
#include <vector>
using namespace std;


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
		//Merge(&M, &S, &j, d_pp[0]*d_pp[2], (d_pp[1]*d_pp[1])*d_pp[2], d_pp[2]);
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



__global__ void MSD_GPU_Interpolate_linear(float *d_MSD_DIT, float *d_MSD_interpolated, int *d_MSD_DIT_widths, int MSD_DIT_size, int *boxcar, int max_width_performed){
	
	int tid  = threadIdx.x;
	if(boxcar[tid] <= max_width_performed) {
	//	int f = threadIdx.x;
		int desired_width = boxcar[tid];
	        int position = (int) floorf(log2f((float) desired_width));
	
	        float width1 = d_MSD_DIT_widths[position];
	        float mean1 = d_MSD_DIT[(position)*MSD_RESULTS_SIZE];
	        float StDev1 = d_MSD_DIT[(position)*MSD_RESULTS_SIZE +1];
	
	//	printf("\nBoxcar: %f \t desired: %f", (float)boxcar[f], desired_width);
	
	        if(position == MSD_DIT_size-1 && width1==(int) desired_width) {
	//                (*mean) = mean1;
	//                (*StDev) = StDev1;
	                  d_MSD_interpolated[tid*2] = mean1;
	                  d_MSD_interpolated[tid*2+1] = StDev1;
	        }
	        else {
	                float width2 = d_MSD_DIT_widths[position+1];
	                float distance_in_width = width2 - width1;
	
	                float mean2 = d_MSD_DIT[(position+1)*MSD_RESULTS_SIZE];
	                float distance_in_mean = mean2 - mean1;
	
	                float StDev2 = d_MSD_DIT[(position+1)*MSD_RESULTS_SIZE +1];
	                float distance_in_StDev = StDev2 - StDev1;
	
	//                        printf("Position: \t %i \t f: %i\n", position, f);
	//                        printf("width:[%f;%f]; mean:[%f;%f]; sd:[%f;%f]\n",width1, width2, mean1, mean2, StDev1, StDev2);
	//                        printf("d width %f; d mean: %f; d StDef: %f\n", distance_in_width, distance_in_mean, distance_in_StDev);
	//                        printf("\tDesired_width: %f\n", desired_width);
	
	//                (*mean) = mean1 + (distance_in_mean/distance_in_width)*((float) desired_width - width1);
	//                (*StDev) = StDev1 + (distance_in_StDev/distance_in_width)*((float) desired_width - width1);
	                d_MSD_interpolated[tid*2] = mean1 + (distance_in_mean/distance_in_width)*((float) desired_width - width1);
	                d_MSD_interpolated[tid*2+1] = StDev1 + (distance_in_StDev/distance_in_width)*((float) desired_width - width1);
	
	        }
	}
}



#endif
