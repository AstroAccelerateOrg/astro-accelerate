#ifndef ASTRO_ACCELERATE_AA_DEVICE_MSD_SHARED_KERNEL_FUNCTIONS_CUH
#define ASTRO_ACCELERATE_AA_DEVICE_MSD_SHARED_KERNEL_FUNCTIONS_CUH

#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include "aa_device_cuda_deprecated_wrappers.cuh"

namespace astroaccelerate {

	//----------------------------------------------------------------------------------------
	//------------- Device functions
	__device__ __inline__ void Initiate(float *M, float *S, float *j, float element) {
		*M = element;
		*S = 0;
		*j = 1.0f;
	}

	__device__ __inline__ void Add_one(float *M, float *S, float *j, float element) {
		float ftemp;
		*j = (*j) + 1.0f;
		*M = (*M) + element;
		ftemp = ((*j)*element - (*M));
		*S = (*S) + 1.0f / ((*j)*((*j) - 1.0f))*ftemp*ftemp;
	}

	__device__ __inline__ void Merge(float *A_M, float *A_S, float *A_j, float B_M, float B_S, float B_j) {
		float ftemp;

		ftemp = (B_j / (*A_j))*(*A_M) - B_M;
		(*A_S) = (*A_S) + B_S + ((*A_j) / (B_j*((*A_j) + B_j)))*ftemp*ftemp;
		(*A_M) = (*A_M) + B_M;
		(*A_j) = (*A_j) + B_j;
	}

	__device__ __inline__ void Reduce_SM(float *M, float *S, float *j, float *s_input) {
		float jv;

		(*M) = s_input[threadIdx.x];
		(*S) = s_input[blockDim.x + threadIdx.x];
		(*j) = s_input[2*blockDim.x + threadIdx.x];

		for (int i = (blockDim.x >> 1); i > HALF_WARP; i = i >> 1) {
			if (threadIdx.x < i) {
				jv = s_input[2*blockDim.x + i + threadIdx.x];
				if (((int)jv)!=0) {
					if ((*j)==0) {
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

	__device__ __inline__ void Reduce_SM_max(float *M, float *S, float *j, float *max, float *min, float *s_par_MSD) {
		int jv;

		(*M) = s_par_MSD[threadIdx.x];
		(*S) = s_par_MSD[blockDim.x + threadIdx.x];
		(*j) = s_par_MSD[2*blockDim.x + threadIdx.x];
		(*max) = s_par_MSD[3*blockDim.x + threadIdx.x];
		(*min) = s_par_MSD[4*blockDim.x + threadIdx.x];


		for (int i = (blockDim.x >> 1); i > HALF_WARP; i = i >> 1) {
			if (threadIdx.x < i) {
				jv = s_par_MSD[2*blockDim.x + i + threadIdx.x];
				if (jv!=0) {
					if ((*j)==0) {
						(*M) = s_par_MSD[i + threadIdx.x];
						(*S) = s_par_MSD[blockDim.x + i + threadIdx.x];
						(*max) = s_par_MSD[3*blockDim.x + i + threadIdx.x];
						(*min) = s_par_MSD[4*blockDim.x + i + threadIdx.x];
						(*j) = jv;
					}
					else {
						Merge(M, S, j, s_par_MSD[i + threadIdx.x], s_par_MSD[blockDim.x + i + threadIdx.x], jv);
						(*max) = fmaxf((*max), s_par_MSD[3*blockDim.x + i + threadIdx.x]);
						(*min) = fminf((*min), s_par_MSD[4*blockDim.x + i + threadIdx.x]);
					}
				}

				s_par_MSD[threadIdx.x] = (*M);
				s_par_MSD[blockDim.x + threadIdx.x] = (*S);
				s_par_MSD[2*blockDim.x + threadIdx.x] = (*j);
				s_par_MSD[3*blockDim.x + threadIdx.x] = (*max);
				s_par_MSD[4*blockDim.x + threadIdx.x] = (*min);
			}
			__syncthreads();
		}
	}

	__device__ __inline__ void Reduce_SM_regular(float *M, float *S, float *j, float *s_input) {
		(*M) = s_input[threadIdx.x];
		(*S) = s_input[blockDim.x + threadIdx.x];
		(*j) = s_input[2*blockDim.x + threadIdx.x];

		for (int i = (blockDim.x >> 1); i > HALF_WARP; i = i >> 1) {
			if (threadIdx.x < i) {
				Merge(M, S, j, s_input[i + threadIdx.x], s_input[blockDim.x + i + threadIdx.x], s_input[2*blockDim.x + i + threadIdx.x]);

				s_input[threadIdx.x] = (*M);
				s_input[blockDim.x + threadIdx.x] = (*S);
				s_input[2*blockDim.x + threadIdx.x] = (*j);
			}
			__syncthreads();
		}
	}

	__device__ __inline__ void Reduce_WARP(float *M, float *S, float *j) {
		float B_M, B_S, B_j;

		for (int q = HALF_WARP; q > 0; q = q >> 1) {
			B_M = aa_shfl_down(AA_ASSUME_MASK, (*M), q);
			B_S = aa_shfl_down(AA_ASSUME_MASK, (*S), q);
			B_j = aa_shfl_down(AA_ASSUME_MASK, (*j), q);

			if (B_j>0) {
				if ((*j)==0) {
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

	__device__ __inline__ void Reduce_WARP_max(float *M, float *S, float *j, float *max, float *min) {
		float B_M, B_S, B_max, B_min, B_j;

		for (int q = HALF_WARP; q > 0; q = q >> 1) {
			B_M = aa_shfl_down(AA_ASSUME_MASK, (*M), q);
			B_S = aa_shfl_down(AA_ASSUME_MASK, (*S), q);
			B_max = aa_shfl_down(AA_ASSUME_MASK, (*max), q);
			B_min = aa_shfl_down(AA_ASSUME_MASK, (*min), q);
			B_j = aa_shfl_down(AA_ASSUME_MASK, (*j), q);

			if (B_j>0) {
				if ((*j)==0) {
					(*S) = B_S;
					(*M) = B_M;
					(*j) = B_j;
					(*max) = B_max;
					(*min) = B_min;
				}
				else {
					Merge(M, S, j, B_M, B_S, B_j);
					(*max) = fmaxf((*max), B_max);
					(*min) = fmaxf((*min), B_min);
				}
			}
		}
	}

	__device__ __inline__ void Reduce_WARP_regular(float *M, float *S, float *j) {
		for (int q = HALF_WARP; q > 0; q = q >> 1) {
			Merge(M, S, j, aa_shfl_down(AA_ASSUME_MASK, (*M), q), aa_shfl_down(AA_ASSUME_MASK, (*S), q), aa_shfl_down(AA_ASSUME_MASK, (*j), q));
		}
	}

	__device__ __inline__ void Sum_partials_regular(float *M, float *S, float *j, float *d_input, float *s_input, int size) {
		int pos;

		//----------------------------------------------
		//---- Summing partials
		pos = threadIdx.x;
		if (size > blockDim.x) {
			(*M) = d_input[MSD_PARTIAL_SIZE*pos];
			(*S) = d_input[MSD_PARTIAL_SIZE*pos + 1];
			(*j) = d_input[MSD_PARTIAL_SIZE*pos + 2];

			pos = pos + blockDim.x;
			while (pos < size) {
				Merge(M, S, j, d_input[MSD_PARTIAL_SIZE*pos], d_input[MSD_PARTIAL_SIZE*pos + 1], d_input[MSD_PARTIAL_SIZE*pos + 2]);
				pos = pos + blockDim.x;
			}

			s_input[threadIdx.x] = (*M);
			s_input[blockDim.x + threadIdx.x] = (*S);
			s_input[2*blockDim.x + threadIdx.x] = (*j);

			__syncthreads();

			Reduce_SM_regular(M, S, j, s_input);
			Reduce_WARP_regular(M, S, j);
		}
		else {
			if (threadIdx.x == 0) {
				pos = 0;
				(*M) = d_input[MSD_PARTIAL_SIZE*pos];
				(*S) = d_input[MSD_PARTIAL_SIZE*pos + 1];
				(*j) = d_input[MSD_PARTIAL_SIZE*pos + 2];

				for (pos = 1; pos < size; pos++) {
					Merge(M, S, j, d_input[MSD_PARTIAL_SIZE*pos], d_input[MSD_PARTIAL_SIZE*pos + 1], d_input[MSD_PARTIAL_SIZE*pos + 2]);
				}
			}
		}
		//---- Summing partials
		//----------------------------------------------
	}

	__device__ __inline__ void Sum_partials_nonregular(float *M, float *S, float *j, float *d_input, float *s_input, int size) {
		int pos;
		float jv;

		//----------------------------------------------
		//---- Summing partials
		pos = threadIdx.x;
		if (size > blockDim.x) {
			(*M) = 0;	(*S) = 0;	(*j) = 0;
			while (pos < size) {
				jv = d_input[MSD_PARTIAL_SIZE*pos + 2];
				if (((int)jv)>0) {
					if ((int)(*j)==0) {
						(*M) = d_input[MSD_PARTIAL_SIZE*pos];
						(*S) = d_input[MSD_PARTIAL_SIZE*pos + 1];
						(*j) = jv;
					}
					else {
						Merge(M, S, j, d_input[MSD_PARTIAL_SIZE*pos], d_input[MSD_PARTIAL_SIZE*pos + 1], jv);
					}
				}
				pos = pos + blockDim.x;
			}

			s_input[threadIdx.x] = (*M);
			s_input[blockDim.x + threadIdx.x] = (*S);
			s_input[2*blockDim.x + threadIdx.x] = (*j);

			__syncthreads();

			Reduce_SM(M, S, j, s_input);
			Reduce_WARP(M, S, j);
		}
		else {
			if (threadIdx.x == 0) {
				pos = 0;
				(*M) = 0;	(*S) = 0;	(*j) = 0;
				for (pos = 1; pos < size; pos++) {
					jv = d_input[MSD_PARTIAL_SIZE*pos + 2];
					if (((int)jv)!=0) {
						if ((int)(*j)==0) {
							(*M) = d_input[MSD_PARTIAL_SIZE*pos];
							(*S) = d_input[MSD_PARTIAL_SIZE*pos + 1];
							(*j) = jv;
						}
						else {
							Merge(M, S, j, d_input[MSD_PARTIAL_SIZE*pos], d_input[MSD_PARTIAL_SIZE*pos + 1], jv);
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
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_MSD_SHARED_KERNEL_FUNCTIONS_CUH
