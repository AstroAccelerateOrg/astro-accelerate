#ifndef DEDISPERSE_KERNEL_H_
#define DEDISPERSE_KERNEL_H_

#define ARRAYSIZE SDIVINT * SDIVINDM

#include "float.h"

// Stores temporary shift values
__device__ __constant__ float dm_shifts[8192];
__device__ __constant__ int   i_nsamp, i_nchans, i_t_processed_s;
__device__ __shared__ ushort2 f_line[UNROLLS][ARRAYSIZE+1];

<<<<<<< HEAD
//{{{ shared_dedisperse_loop

//__launch_bounds__(SDIVINT*SDIVINDM)
=======
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
__global__ void shared_dedisperse_kernel(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	int	i, j, c;
	int	shift[UNROLLS];
	
	ushort	temp_f;
	int	local, unroll;
	
<<<<<<< HEAD
	float 	findex = (threadIdx.x*2);
	float	local_kernel_one[SNUMREG];
	float	local_kernel_two[SNUMREG];
=======
	float findex = (threadIdx.x*2);
	float local_kernel_one[SNUMREG];
	float local_kernel_two[SNUMREG];
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da

	for(i = 0; i < SNUMREG; i++) {
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int	idx = (threadIdx.x + (threadIdx.y * SDIVINT));
	int	nsamp_counter = (idx + (blockIdx.x*(2*SNUMREG*SDIVINT)));

	float	shift_two = (mstartdm + (__int2float_rz(blockIdx.y)*SFDIVINDM*mdmstep));
	float	shift_one = (__int2float_rz(threadIdx.y)*mdmstep);

<<<<<<< HEAD
	for(c = 0; c < i_nchans; c+=UNROLLS) {

		__syncthreads();
	
		for(j=0; j<UNROLLS; j++) { 
			temp_f = (__ldg((d_input +( __float2int_rz(dm_shifts[c+j]*shift_two)))+ (nsamp_counter + (j*i_nsamp))));

			f_line[j][idx].x = temp_f;
			if(idx > 0) {
=======
	for(c = 0; c < i_nchans; c+=UNROLLS)
	{
		__syncthreads();
	
		for(j=0; j<UNROLLS; j++) 
		{ 
			temp_f = (__ldg((d_input +( __float2int_rz(dm_shifts[c+j]*shift_two)))+ (nsamp_counter + (j*i_nsamp))));

			f_line[j][idx].x = temp_f;
			if(idx > 0) 
			{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
				f_line[j][idx-1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one*dm_shifts[c+j] + findex);
		}

		nsamp_counter = (nsamp_counter + (UNROLLS*i_nsamp));

		__syncthreads();

<<<<<<< HEAD
		for(i = 0; i < SNUMREG; i++) {
			local = 0;
			unroll = (i*2*SDIVINT);
			for(j=0; j<UNROLLS; j++) {
=======
		for(i = 0; i < SNUMREG; i++) 
		{
			local = 0;
			unroll = (i*2*SDIVINT);
			for(j=0; j<UNROLLS; j++) 
			{
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
				local += *(int*)&f_line[j][(shift[j] + unroll)];
			}
			local_kernel_one[i] += ((ushort2*)(&local))->x;
			local_kernel_two[i] += ((ushort2*)(&local))->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ((((blockIdx.y*SDIVINDM) + threadIdx.y)*(i_t_processed_s)) + (blockIdx.x*2*SNUMREG*SDIVINT)) + 2*threadIdx.x;

	#pragma unroll
<<<<<<< HEAD
	for(i = 0; i < SNUMREG; i++) {
		*((float2*)(d_output+local + (i*2*SDIVINT))) = make_float2(local_kernel_one[i]/i_nchans/bin, local_kernel_two[i]/i_nchans/bin);
//		d_output[local + (i*2*SDIVINT)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT) + 1] = (local_kernel_two[i])/i_nchans;
=======
	for(i = 0; i < SNUMREG; i++) 
	{
		*((float2*)(d_output+local + (i*2*SDIVINT))) = make_float2(local_kernel_one[i]/i_nchans/bin, local_kernel_two[i]/i_nchans/bin);
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	}
}

#endif

