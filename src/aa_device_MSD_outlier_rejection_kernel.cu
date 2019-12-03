#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "aa_params.hpp"
#include "aa_device_MSD_outlier_rejection_kernel.hpp"
#include "aa_device_MSD_shared_kernel_functions.cuh"
#include "aa_device_cuda_deprecated_wrappers.cuh"

namespace astroaccelerate {

	__global__ void MSD_GPU_calculate_partials_2d_and_minmax(float const* __restrict__ d_input, float *d_output, int y_steps, int nTimesamples, int offset) {
		__shared__ float s_input[5*PD_NTHREADS];
		float M, S, max, min, ftemp, j;

		int spos = blockIdx.x*PD_NTHREADS + threadIdx.x;
		size_t gpos = blockIdx.y*(size_t)y_steps*(size_t)nTimesamples + (size_t)spos;
		M=0;	S=0;	j=0;	max=0;	min=0;
		if( spos<(nTimesamples-offset) ){	
			ftemp = (float) d_input[gpos];
			Initiate( &M, &S, &j, ftemp);
			max = ftemp;
			min = ftemp;
			
			gpos = gpos + (size_t)nTimesamples;
			for (int yf = 1; yf < y_steps; yf++) {
				ftemp = (float) d_input[gpos];
				max = (fmaxf(max,ftemp));
				min = (fminf(min,ftemp));
				Add_one( &M, &S, &j, ftemp);
				gpos = gpos + (size_t)nTimesamples;
			}
		}
		
		s_input[threadIdx.x] = M;
		s_input[blockDim.x + threadIdx.x] = S;
		s_input[2*blockDim.x + threadIdx.x] = j;
		s_input[3*blockDim.x + threadIdx.x] = max;
		s_input[4*blockDim.x + threadIdx.x] = min;
		
		__syncthreads();
		
		Reduce_SM_max( &M, &S, &j, &max, &min, s_input );
		Reduce_WARP_max( &M, &S, &j, &max, &min);
		
		//----------------------------------------------
		//---- Writing data
		if (threadIdx.x == 0) {
			gpos = blockIdx.y*gridDim.x + blockIdx.x;
			d_output[MSD_PARTIAL_SIZE*gpos] = M;
			d_output[MSD_PARTIAL_SIZE*gpos + 1] = S;
			d_output[MSD_PARTIAL_SIZE*gpos + 2] = j;
			d_output[MSD_PARTIAL_SIZE*gpos + 3] = max;
			d_output[MSD_PARTIAL_SIZE*gpos + 4] = min;
			//if(blockIdx.x<10 && blockIdx.y<10) printf("Initial: bl:[%d;%d]; M=%f; S=%f; j=%f; max=%f; min=%f\n", blockIdx.x, blockIdx.y, M, S, j, max, min);
		}
	}

	__global__ void MSD_BLN_calculate_partials_2d_and_minmax_with_outlier_rejection(float const* __restrict__ d_input, float *d_output, float *d_MSD, int y_steps, int nTimesamples, int offset, float bln_sigma_constant) {
		__shared__ float s_input[5*PD_NTHREADS];
		float M, S, j, ftemp, max, min;
		float limit_down = d_MSD[0] - bln_sigma_constant*d_MSD[1];
		float limit_up   = d_MSD[0] + bln_sigma_constant*d_MSD[1];

		size_t temp_gpos = blockIdx.y*gridDim.x + blockIdx.x;
		max = d_output[MSD_PARTIAL_SIZE*temp_gpos + 3];
		min = d_output[MSD_PARTIAL_SIZE*temp_gpos + 4];
		if( (min>limit_down) && (max < limit_up) ) return;
	
		int spos = blockIdx.x*PD_NTHREADS + threadIdx.x;
		size_t gpos = blockIdx.y*(size_t)y_steps*(size_t)nTimesamples + (size_t)spos;
		M=0;	S=0;	j=0;	max=0;	min=0;
		if( spos<(nTimesamples-offset) ){
			for (int yf = 0; yf < y_steps; yf++) {
				ftemp = (float) d_input[gpos];
				if( (ftemp>limit_down) && (ftemp < limit_up) ){
					if(j==0){
						Initiate( &M, &S, &j, ftemp);
						max = ftemp;
						min = ftemp;
					}
					else{
						Add_one( &M, &S, &j, ftemp);
						max = fmaxf(max, ftemp);
						min = fminf(min, ftemp);
					}			
				}
				gpos = gpos + (size_t)nTimesamples;
			}
		}

		s_input[threadIdx.x] = M;
		s_input[blockDim.x + threadIdx.x] = S;
		s_input[2*blockDim.x + threadIdx.x] = j;
		s_input[3*blockDim.x + threadIdx.x] = max;
		s_input[4*blockDim.x + threadIdx.x] = min;

		__syncthreads();

		Reduce_SM_max( &M, &S, &j, &max, &min, s_input );
		Reduce_WARP_max( &M, &S, &j, &max, &min);

		//----------------------------------------------
		//---- Writing data
		if (threadIdx.x == 0) {
			gpos = blockIdx.y*gridDim.x + blockIdx.x;
			d_output[MSD_PARTIAL_SIZE*gpos] = M;
			d_output[MSD_PARTIAL_SIZE*gpos + 1] = S;
			d_output[MSD_PARTIAL_SIZE*gpos + 2] = j;
			d_output[MSD_PARTIAL_SIZE*gpos + 3] = max;
			d_output[MSD_PARTIAL_SIZE*gpos + 4] = min;
			//if(blockIdx.x<10 && blockIdx.y<10) printf("bl:[%d;%d]; M=%f; S=%f; j=%f; max=%f; min=%f\n", blockIdx.x, blockIdx.y, M, S, j, max, min);
		}
	}


	__global__ void MSD_BLN_grid_outlier_rejection_new(float *d_input, float *d_output, int size, float multiplier) {
		__shared__ float s_input[3*WARP*WARP];
		__shared__ float s_signal_mean;
		__shared__ float s_signal_sd;
		float M, S, j;
		float signal_mean, signal_sd;

		//----------------------------------------------
		//---- Calculation of the initial MSD
		Sum_partials_nonregular(&M, &S, &j, d_input, s_input, size);

		if (threadIdx.x==0) {
			s_signal_mean = M/j;
			s_signal_sd = sqrt(S/j);
		}

		__syncthreads();

		signal_mean = s_signal_mean;
		signal_sd = s_signal_sd;
		//---- Calculation of the initial MSD
		//----------------------------------------------


		//----------------------------------------------
		//---- Iterations with outlier rejection
		for (int f = 0; f<5; f++) {
			int pos;
			float jv, Mt;

			pos = threadIdx.x;
			if (size > blockDim.x) {
				M = 0;	S = 0;	j = 0;
				while (pos < size) {
					jv = __ldg(&d_input[3*pos + 2]);
					Mt = __ldg(&d_input[3*pos]);
					if (((int)jv)!=0) {
						if ((Mt/jv >(signal_mean - multiplier*signal_sd)) && (Mt/jv < (signal_mean + multiplier*signal_sd))) {
							if ((int)j==0) {
								M = Mt;
								S = __ldg(&d_input[3*pos + 1]);
								j = jv;
							}
							else {
								Merge(&M, &S, &j, Mt, __ldg(&d_input[3*pos + 1]), jv);
							}
						}
					}
					pos = pos + blockDim.x;
				}

				s_input[threadIdx.x] = M;
				s_input[blockDim.x + threadIdx.x] = S;
				s_input[2*blockDim.x + threadIdx.x] = j;

				__syncthreads();

				Reduce_SM(&M, &S, &j, s_input);
				Reduce_WARP(&M, &S, &j);
			}
			else {
				if (threadIdx.x == 0) {
					pos = 0;
					M = 0; S = 0; j = 0;
					for (pos = 0; pos < size; pos++) {
						jv = __ldg(&d_input[3*pos + 2]);
						Mt = __ldg(&d_input[3*pos]);
						if (((int)jv)!=0) {
							if ((Mt/jv >(signal_mean - multiplier*signal_sd)) && (Mt/jv < (signal_mean + multiplier*signal_sd))) {
								if ((int)j==0) {
									M = Mt;
									S = __ldg(&d_input[3*pos + 1]);
									j = jv;
								}
								else {
									Merge(&M, &S, &j, Mt, __ldg(&d_input[3*pos + 1]), jv);
								}
							}
						}
					}
				}
			}

			if (threadIdx.x == 0) {
				s_signal_mean = M/j;
				s_signal_sd = sqrt(S/j);
			}

			__syncthreads();

			signal_mean = s_signal_mean;
			signal_sd = s_signal_sd;
		}
		//---- Iterations with outlier rejection
		//----------------------------------------------



		//----------------------------------------------
		//---- Writing data
		if (threadIdx.x==0) {
			d_output[0] = signal_mean;
			d_output[1] = signal_sd;
			d_output[2] = j;
		}
		//---- Writing data
		//----------------------------------------------
	}


	__global__ void MSD_BLN_grid_calculate_partials(float const* __restrict__ d_input, float *d_output, int x_steps, int y_steps, int nColumns, int msd) {
		extern __shared__ float Ms_Ss[];

		int warp_id, local_id, dim_y;
		size_t pos;
		float x; // current element
		float M; // streaming mean
		float S; // streaming sum of squares (x_i-\bar{x})
		float j;
		float ftemp;

		local_id = threadIdx.x & (WARP - 1);
		warp_id = threadIdx.x>>5;
		dim_y = blockDim.x>>5;

		//----------------------------------------------
		//---- Calculating of streaming mean and sum of squares
		pos = (blockIdx.y*(size_t)dim_y + (size_t)warp_id)*(size_t)y_steps*(size_t)nColumns + blockIdx.x*WARP*(size_t)x_steps + (size_t)local_id;
		M = __ldg(&d_input[pos]);
		S = 0;
		j = 1.0f;
		for (int xf = 1; xf<x_steps; xf++) {
			pos = pos + WARP;
			x = __ldg(&d_input[pos]);
			j = j+1.0f;
			M = M + x;
			ftemp = (j*x - M);
			S = S + 1.0f/(j*(j-1.0f))*ftemp*ftemp;
		}

		pos = pos + (size_t)nColumns - (size_t)(x_steps-1)*WARP;
		for (int yf = 1; yf<y_steps; yf++) {
			for (int xf = 0; xf<x_steps; xf++) {
				x = __ldg(&d_input[pos]);
				j = j+1.0f;
				M = M + x;
				ftemp = (j*x - M);
				S = S + 1.0f/(j*(j-1.0f))*ftemp*ftemp;
				pos = pos + WARP;
			}
			pos = pos + (size_t)nColumns - (size_t)x_steps*WARP;
		}

		Ms_Ss[threadIdx.x] = M;
		Ms_Ss[blockDim.x + threadIdx.x] = S;

		__syncthreads();

		// now all threads had saved their work, reduction follows

		// first we must load initial values
		//j=Neco;
		for (int i = (blockDim.x>>1); i>HALF_WARP; i = i>>1) {
			if (threadIdx.x<i) {
				j = j*2;
				ftemp = (M - Ms_Ss[i + threadIdx.x]);
				S = S + Ms_Ss[blockDim.x + i + threadIdx.x] + (1.0f/j)*ftemp*ftemp;
				M = M + Ms_Ss[i + threadIdx.x];

				Ms_Ss[threadIdx.x] = M;
				Ms_Ss[blockDim.x + threadIdx.x] = S;
			}
			// in the last iteration we do not need to save the results... or synchronize threads...
			__syncthreads();
		}

		// by now we should have only 32 partial results. shuffle reduction follows
		for (int q = HALF_WARP; q>0; q = q>>1) {
			j = j*2;
			ftemp = (M - aa_shfl_down(AA_ASSUME_MASK, M, q));
			S = S + aa_shfl_down(AA_ASSUME_MASK, S, q) + (1.0f/j)*ftemp*ftemp;
			M = M + aa_shfl_down(AA_ASSUME_MASK, M, q);
		}

		//----------------------------------------------
		//---- Writing data
		if (threadIdx.x==0) {
			pos = blockIdx.y*gridDim.x + blockIdx.x;
			if (msd) {
				// produce mean and sd instead of T and S
				d_output[3*pos] = M/j;
				d_output[3*pos + 1] = sqrt(S/j);
			}
			else {
				d_output[3*pos] = M;
				d_output[3*pos + 1] = S;
			}
		}
	}


	__global__ void MSD_BLN_grid_outlier_rejection(float const* __restrict__ d_input, float *d_output, int size, float nElements, float multiplier) {
		__shared__ float Ms[WARP*WARP];
		__shared__ float Ss[WARP*WARP];
		__shared__ float js[WARP*WARP];
		__shared__ float s_signal_mean;
		__shared__ float s_signal_sd;


		int  pos; //warp_id,
		float M, Mt, S, j, jv;
		float ftemp;
		float signal_mean, signal_sd;

		//warp_id = threadIdx.x>>5;

		//----------------------------------------------
		//---- Calculation of the initial MSD
		pos = threadIdx.x;
		if (size>blockDim.x) {
			M = __ldg(&d_input[3*pos]);
			S = __ldg(&d_input[3*pos+1]);
			j = nElements;
			pos = pos + blockDim.x;
			while (pos<size) {
				jv = nElements;
				ftemp = (jv/j*M - __ldg(&d_input[3*pos]));
				S = S + __ldg(&d_input[3*pos+1]) + (j/(jv*(j+jv)))*ftemp*ftemp;
				M = M + __ldg(&d_input[3*pos]);
				j = j+jv;
				pos = pos + blockDim.x;
			}

			__syncthreads();

			Ms[threadIdx.x] = M;
			Ss[threadIdx.x] = S;
			js[threadIdx.x] = j;
			// now all threads had saved their work, reduction follows		
			// first we must load initial values
			for (int i = (blockDim.x>>1); i>HALF_WARP; i = i>>1) {
				if (threadIdx.x<i) {
					jv = js[i + threadIdx.x];
					ftemp = (jv/j*M - Ms[i + threadIdx.x]);
					S = S + Ss[i + threadIdx.x] + (j/(jv*(j+jv)))*ftemp*ftemp;
					M = M + Ms[i + threadIdx.x];
					j = j+jv;

					Ms[threadIdx.x] = M;
					Ss[threadIdx.x] = S;
					js[threadIdx.x] = j;
				}
				__syncthreads();
			}

			// by now we should have only 32 partial results. shuffle reduction follows
			for (int q = HALF_WARP; q>0; q = q>>1) {
				jv = aa_shfl_down(AA_ASSUME_MASK, j, q);
				ftemp = (jv/j*M - aa_shfl_down(AA_ASSUME_MASK, M, q));
				S = S + aa_shfl_down(AA_ASSUME_MASK, S, q) + (j/(jv*(j+jv)))*ftemp*ftemp;
				M = M + aa_shfl_down(AA_ASSUME_MASK, M, q);
				j = j+jv;
			}

		}
		else {
			if (threadIdx.x==0) {
				pos = 0;
				M = __ldg(&d_input[3*pos]);
				S = __ldg(&d_input[3*pos+1]);
				j = nElements;
				for (pos = 1; pos<size; pos++) {
					jv = __ldg(&d_input[3*pos+2]);
					ftemp = (jv/j*M - __ldg(&d_input[3*pos]));
					S = S + __ldg(&d_input[3*pos+1]) + (j/(jv*(j+jv)))*ftemp*ftemp;
					M = M + __ldg(&d_input[3*pos]);
					j = j+jv;
				}
			}
		}

		if (threadIdx.x==0) {
			s_signal_mean = M/j;
			s_signal_sd = sqrt(S/j);
		}

		__syncthreads();

		signal_mean = s_signal_mean;
		signal_sd = s_signal_sd;
		//---- Calculation of the initial MSD
		//----------------------------------------------

		//if(threadIdx.x==0) printf("Initial mean:%f; and standard deviation:%f;\n", signal_mean, signal_sd);

		//----------------------------------------------
		//---- Iterations with outlier rejection
		for (int f = 0; f<5; f++) {
			pos = threadIdx.x;
			if (size>blockDim.x) {
				M = 0;
				S = 0;
				j = 0;
				while (pos<size) {
					Mt = __ldg(&d_input[3*pos]);
					if ((Mt/nElements >(signal_mean - multiplier*signal_sd)) && (Mt/nElements < (signal_mean + multiplier*signal_sd))) {
						if (j==0) {
							M = Mt;
							S = __ldg(&d_input[3*pos+1]);
							j = nElements;
						}
						else {
							jv = nElements;
							ftemp = (jv/j*M - Mt);
							S = S + __ldg(&d_input[3*pos+1]) + (j/(jv*(j+jv)))*ftemp*ftemp;
							M = M + Mt;
							j = j+jv;
						}
					}
					pos = pos + blockDim.x;
				}

				__syncthreads();

				Ms[threadIdx.x] = M;
				Ss[threadIdx.x] = S;
				js[threadIdx.x] = j;
				// now all threads had saved their work, reduction follows		
				// first we must load initial values
				for (int i = (blockDim.x>>1); i>HALF_WARP; i = i>>1) {
					if (threadIdx.x<i) {
						jv = js[i + threadIdx.x];
						if (jv!=0) {
							if (j==0) {
								S = Ss[i + threadIdx.x];
								M = Ms[i + threadIdx.x];
								j = jv;
							}
							else {
								ftemp = (jv/j*M - Ms[i + threadIdx.x]);
								S = S + Ss[i + threadIdx.x] + (j/(jv*(j+jv)))*ftemp*ftemp;
								M = M + Ms[i + threadIdx.x];
								j = j+jv;
							}
						}

						Ms[threadIdx.x] = M;
						Ss[threadIdx.x] = S;
						js[threadIdx.x] = j;
					}
					__syncthreads();
				}

				// by now we should have only 32 partial results. shuffle reduction follows
				for (int q = HALF_WARP; q>0; q = q>>1) {
					jv = aa_shfl_down(AA_ASSUME_MASK, j, q);
					if (jv!=0) {
						if (j==0) {
							S = aa_shfl_down(AA_ASSUME_MASK, S, q);
							M = aa_shfl_down(AA_ASSUME_MASK, M, q);
							j = jv;
						}
						else {
							ftemp = (jv/j*M - aa_shfl_down(AA_ASSUME_MASK, M, q));
							S = S + aa_shfl_down(AA_ASSUME_MASK, S, q) + (j/(jv*(j+jv)))*ftemp*ftemp;
							M = M + aa_shfl_down(AA_ASSUME_MASK, M, q);
							j = j+jv;
						}

					}
				}

			}
			else {
				if (threadIdx.x==0) {
					M = 0;
					S = 0;
					j = 0;
					for (pos = 0; pos<size; pos++) {
						Mt = __ldg(&d_input[3*pos]);
						if ((Mt/nElements >(signal_mean - multiplier*signal_sd)) && (Mt/nElements < (signal_mean + multiplier*signal_sd))) {
							if (j==0) {
								M = Mt;
								S = __ldg(&d_input[3*pos+1]);
								j = nElements;
							}
							else {
								jv = nElements;
								ftemp = (jv/j*M - __ldg(&d_input[3*pos]));
								S = S + __ldg(&d_input[3*pos+1]) + (j/(jv*(j+jv)))*ftemp*ftemp;
								M = M + __ldg(&d_input[3*pos]);
								j = j+jv;
							}
						}
					}
				}
			}

			if (threadIdx.x==0) {
				s_signal_mean = M/j;
				s_signal_sd = sqrt(S/j);
			}

			__syncthreads();

			signal_mean = s_signal_mean;
			signal_sd = s_signal_sd;

			//if(threadIdx.x==0) printf("Corrected mean:%f; and standard deviation:%f;\n", signal_mean, signal_sd);
		}
		//---- Iterations with outlier rejection
		//----------------------------------------------



		//----------------------------------------------
		//---- Writing data
		if (threadIdx.x==0) {
			d_output[0] = signal_mean;
			d_output[1] = signal_sd;
			d_output[2] = j;
		}
		//---- Writing data
		//----------------------------------------------
	}



	//----------------------------------------------------------------
	//----------------------------------------------------------------
	//----------------------------------------------------------------

	/** \brief Wrapper function to MSD_BLN_pw_rejection_normal kernel function. */
	void call_kernel_MSD_BLN_calculate_partials_2d_and_minmax_with_outlier_rejection(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, float *const d_MSD, const int &y_steps, const int &nTimesamples, const int &offset, const float &bln_sigma_constant) {
		MSD_BLN_calculate_partials_2d_and_minmax_with_outlier_rejection<<<grid_size, block_size>>>(d_input, d_output, d_MSD, y_steps, nTimesamples, offset, bln_sigma_constant);
	}
	
	void call_kernel_MSD_GPU_calculate_partials_2d_and_minmax(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, const int &y_steps, const int &nTimesamples, const int &offset) {
		MSD_GPU_calculate_partials_2d_and_minmax<<<grid_size, block_size>>>(d_input, d_output, y_steps, nTimesamples, offset);
	}

	/** \brief Wrapper function to MSD_BLN_grid_outlier_rejection_new kernel function. */
	void call_kernel_MSD_BLN_grid_outlier_rejection_new(const dim3 &grid_size, const dim3 &block_size, float *const d_input, float *const d_output, const int &size, const float &multiplier) {
		MSD_BLN_grid_outlier_rejection_new<<<grid_size, block_size>>>(d_input, d_output, size, multiplier);
	}

	/** \brief Wrapper function to MSD_BLN_grid_calculate_partials kernel function. */
	void call_kernel_MSD_BLN_grid_calculate_partials(const dim3 &grid_size, const dim3 &block_size, const int &threads, float const *const d_input, float *const d_output, const int &x_steps, const int &y_steps, const int &nColumns, const int &msd) {
		MSD_BLN_grid_calculate_partials<<<grid_size, block_size, threads>>>(d_input, d_output, x_steps, y_steps, nColumns, msd);
	}

	/** \brief Wrapper function to MSD_BLN_grid_outlier_rejection kernel function. */
	void call_kernel_MSD_BLN_grid_outlier_rejection(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, const int &size, const float &nElements, const float &multiplier) {
		MSD_BLN_grid_outlier_rejection<<<grid_size, block_size>>>(d_input, d_output, size, nElements, multiplier);
	}

} //namespace astroaccelerate
