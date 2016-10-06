#ifndef DEDISPERSE_KERNEL_H_
#define DEDISPERSE_KERNEL_H_

#define ARRAYSIZE SDIVINT * SDIVINDM

<<<<<<< HEAD
#include "float.h"

// Stores temporary shift values
__device__ __constant__ float dm_shifts[8192];
__device__ __constant__ int   i_nsamp, i_nchans, i_t_processed_s;
__device__ __shared__ ushort2 f_line[UNROLLS][ARRAYSIZE+1];

__global__ void shared_dedisperse_kernel(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	int	i, j, c;
	int	shift[UNROLLS];
	
	ushort	temp_f;
	int	local, unroll;
	
	float findex = (threadIdx.x*2);
	float local_kernel_one[SNUMREG];
	float local_kernel_two[SNUMREG];

	for(i = 0; i < SNUMREG; i++) {
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int	idx = (threadIdx.x + (threadIdx.y * SDIVINT));
	int	nsamp_counter = (idx + (blockIdx.x*(2*SNUMREG*SDIVINT)));

	float	shift_two = (mstartdm + (__int2float_rz(blockIdx.y)*SFDIVINDM*mdmstep));
	float	shift_one = (__int2float_rz(threadIdx.y)*mdmstep);

	for(c = 0; c < i_nchans; c+=UNROLLS)
	{
		__syncthreads();
	
		for(j=0; j<UNROLLS; j++) 
		{ 
			temp_f = (__ldg((d_input +( __float2int_rz(dm_shifts[c+j]*shift_two)))+ (nsamp_counter + (j*i_nsamp))));

			f_line[j][idx].x = temp_f;
			if(idx > 0) 
			{
				f_line[j][idx-1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one*dm_shifts[c+j] + findex);
		}

		nsamp_counter = (nsamp_counter + (UNROLLS*i_nsamp));

		__syncthreads();

		for(i = 0; i < SNUMREG; i++) 
		{
			local = 0;
			unroll = (i*2*SDIVINT);
			for(j=0; j<UNROLLS; j++) 
			{
				local += *(int*)&f_line[j][(shift[j] + unroll)];
			}
			local_kernel_one[i] += ((ushort2*)(&local))->x;
			local_kernel_two[i] += ((ushort2*)(&local))->y;
=======
// Stores temporary shift values
__device__ __constant__ float dm_shifts[15500];
__device__ __constant__ int   i_nsamp, i_nchans, i_t_processed_s;
__device__ __shared__ float2 fa_line[ARRAYSIZE+1];
__device__ __shared__ float2 fb_line[ARRAYSIZE+1];
__device__ __shared__ float2 fc_line[ARRAYSIZE+1];
__device__ __shared__ float2 fd_line[ARRAYSIZE+1];

//{{{ shared_dedisperse_loop

__global__ void shared_dedisperse_kernel(float *d_input, float *d_output, cudaTextureObject_t tex, float mstartdm, float mdmstep)
{
	// NOTE: inshift AND outshift are set to 0 (zero) in the kernel call and so is
	// removed from this kernel.
	
	int   shift;	
	float local_kernel_t[SNUMREG];

	int t  = blockIdx.x * SNUMREG * SDIVINT  + threadIdx.x;
	
	// Initialise the time accumulators
	for(int i = 0; i < SNUMREG; i++) local_kernel_t[i] = 0.0f;

	float shift_temp = mstartdm + ((blockIdx.y * SDIVINDM + threadIdx.y) * mdmstep);
	
	// Loop over the frequency channels.
        for(int c = 0; c < i_nchans; c++) {


		// Calculate the initial shift for this given frequency
		// channel (c) at the current despersion measure (dm) 
		// ** dm is constant for this thread!!**
		shift = (c * (i_nsamp) + t) + __float2int_rz (dm_shifts[c] * shift_temp);
		
		#pragma unroll
		for(int i = 0; i < SNUMREG; i++) {
#ifdef SM_35
			local_kernel_t[i] += __ldg(d_input + shift + (i * SDIVINT));
#else
			local_kernel_t[i] += d_input[shift + (i * SDIVINT) ];
#endif			
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
		}
	}

	// Write the accumulators to the output array. 
<<<<<<< HEAD
	local = ((((blockIdx.y*SDIVINDM) + threadIdx.y)*(i_t_processed_s)) + (blockIdx.x*2*SNUMREG*SDIVINT)) + 2*threadIdx.x;

	#pragma unroll
	for(i = 0; i < SNUMREG; i++) 
	{
		*((float2*)(d_output+local + (i*2*SDIVINT))) = make_float2(local_kernel_one[i]/i_nchans/bin, local_kernel_two[i]/i_nchans/bin);
	}
}

=======
	#pragma unroll
	for(int i = 0; i < SNUMREG; i++) {
		d_output[((blockIdx.y * SDIVINDM) + threadIdx.y)*(i_t_processed_s) + (i * SDIVINT) + (SNUMREG * SDIVINT * blockIdx.x) + threadIdx.x] = local_kernel_t[i]/i_nchans;
	}
}
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
#endif

