#ifndef DEDISPERSE_KERNEL_H_
#define DEDISPERSE_KERNEL_H_

#define ARRAYSIZE DIVINT * DIVINDM

#include <cuda.h>
#include <cuda_runtime.h>
#include <cutil_inline.h>
#include "params.h"

// Stores temporary shift values
//__device__ __constant__ float dm_shifts[8192];
__device__ __constant__ float dm_shifts[16000];
__device__ __constant__ int   i_nsamp, i_maxshift, i_nchans;
__device__ __shared__ float f_line[ARRAYSIZE];

//{{{ cache_dedisperse_loop
__global__ void cache_dedisperse_loop(float *outbuff, float *buff, float mstartdm, float mdmstep)
{

	// NOTE: inshift AND outshift are set to 0 (zero) in the kernel call and so is
	// removed from this kernel.
	
	int   shift;	
	float local_kernel_t[NUMREG];

	int t  = blockIdx.x * NUMREG * DIVINT  + threadIdx.x;
	
	// Initialise the time accumulators
	for(int i = 0; i < NUMREG; i++) local_kernel_t[i] = 0.0f;

	float shift_temp = mstartdm + ((blockIdx.y * DIVINDM + threadIdx.y) * mdmstep);
	
	// Loop over the frequency channels.
        for(int c = 0; c < i_nchans; c++) {


		// Calculate the initial shift for this given frequency
		// channel (c) at the current despersion measure (dm) 
		// ** dm is constant for this thread!!**
		shift = (c * (i_nsamp) + t) + __float2int_rz (dm_shifts[c] * shift_temp);
		
		#pragma unroll
		for(int i = 0; i < NUMREG; i++) {
			local_kernel_t[i] += buff[shift + (i * DIVINT) ];
		}
	}

	// Write the accumulators to the output array. 
	#pragma unroll
	for(int i = 0; i < NUMREG; i++) {
		outbuff[((blockIdx.y * DIVINDM) + threadIdx.y)* (i_nsamp-i_maxshift) + (i * DIVINT) + (NUMREG * DIVINT * blockIdx.x) + threadIdx.x] = local_kernel_t[i];
	}

}

//}}}

//{{{ cache_contiguous_loop
__global__ void cache_contiguous_loop(float *outbuff, float *buff, float mstartdm, float mdmstep)
{

	// NOTE: inshift AND outshift are set to 0 (zero) in the kernel call and so is
	// removed from this kernel.
	
	int   shift;	
	float local_kernel_t[NUMREG];

	//int t  = blockIdx.x * NUMREG * DIVINT  + threadIdx.x;
	int t  = blockIdx.x * DIVINT * NUMREG + NUMREG * threadIdx.x;
	
	// Initialise the time accumulators
	for(int i = 0; i < NUMREG; i++) local_kernel_t[i] = 0.0f;

	float shift_temp = mstartdm + ((blockIdx.y * DIVINDM + threadIdx.y) * mdmstep);
	
	// Loop over the frequency channels.
        for(int c = 0; c < i_nchans; c++) {


		// Calculate the initial shift for this given frequency
		// channel (c) at the current despersion measure (dm) 
		// ** dm is constant for this thread!!**
		shift = (c * (i_nsamp) + t) + __float2int_rz (dm_shifts[c] * shift_temp);
		
		#pragma unroll
		for(int i = 0; i < NUMREG; i++) {
			local_kernel_t[i] += buff[shift + i];
		}
	}

	// Write the accumulators to the output array. 
	#pragma unroll
	for(int i = 0; i < NUMREG; i++) {
		outbuff[((blockIdx.y * DIVINDM) + threadIdx.y)* (i_nsamp-i_maxshift) + (i) + t] = local_kernel_t[i];
	}

}

//}}}

//{{{ shared_dedisperse_loop

__global__ void shared_dedisperse_loop(float *outbuff, float *buff, float mstartdm, float mdmstep)
{
	int   i, c, shift;	
	float local_kernel_t[NUMREG];

	float by = blockIdx.y*DIVINDM*mdmstep;
	float ty = threadIdx.y;
	
	float bx = blockIdx.x*NUMREG*DIVINT;
	float tx = threadIdx.x;
	
	int idx = (threadIdx.x + (threadIdx.y * DIVINT));

	float shift_one = (mstartdm +(by + (ty*mdmstep)));
	float shift_two = (mstartdm + by);

	// Initialise the time accumulators
	
	#pragma unroll
	for(i = 0; i < NUMREG; i++) local_kernel_t[i] = 0.0f;
        
	for(c = 0; c < i_nchans; c++) {

		//f_line[idx] = buff[((c*i_nsamp) + idx) + __float2int_rz(((dm_shifts[c]*shift_one) + bx))];
		f_line[idx] = buff[((c*i_nsamp) + idx) + __float2int_rz(((dm_shifts[c]*shift_two) + bx))];
		shift = __float2int_rz(((dm_shifts[c]*shift_one) + tx) - (dm_shifts[c]*shift_two));
		__syncthreads();

		#pragma unroll
		for(i = 0; i < NUMREG; i++) local_kernel_t[i] += f_line[(shift + (i*DIVINT))];

	}

	// Write the accumulators to the output array. 
	shift = ((blockIdx.y*DIVINDM) + threadIdx.y)*(i_nsamp-i_maxshift) + (int)bx + threadIdx.x;

	#pragma unroll
	for(i = 0; i < NUMREG; i++) {
		outbuff[shift + (i*DIVINT)] = local_kernel_t[i];
	}
}

//}}}

//{{{ shared_contiguous_loop

__global__ void shared_contiguous_loop(float *outbuff, float *buff, float mstartdm, float mdmstep)
{

	int   i, c, shift;	
	float local_kernel_t[NUMREG];

	// Initialise the time accumulators
	for(i = 0; i < NUMREG; i++) local_kernel_t[i] = 0.0f;

	int shift_one = (mstartdm +((blockIdx.y*DIVINDM + threadIdx.y)*mdmstep));
	int shift_two = (mstartdm + (blockIdx.y*DIVINDM*mdmstep));
	int idx = (threadIdx.x + (threadIdx.y * DIVINT));
        for(c = 0; c < i_nchans; c++) {
		
		f_line[idx] = buff[((c*i_nsamp) + (blockIdx.x*NUMREG*DIVINT + idx)) + __float2int_rz(dm_shifts[c]*shift_two)];
		__syncthreads();
		
		shift = __float2int_rz(dm_shifts[c]*shift_one) - __float2int_rz(dm_shifts[c]*shift_two) + threadIdx.x*NUMREG;

		#pragma unroll
		for(i = 0; i < NUMREG; i++) {
			local_kernel_t[i] += f_line[(shift + i)];
		}
	}

	// Write the accumulators to the output array. 
	shift = ((blockIdx.y*DIVINDM) + threadIdx.y)*(i_nsamp-i_maxshift) + (NUMREG*DIVINT*blockIdx.x) + threadIdx.x*NUMREG;

	#pragma unroll
	for(i = 0; i < NUMREG; i++) {
		outbuff[shift + i] = local_kernel_t[i];
	}
}

//}}}

#endif
