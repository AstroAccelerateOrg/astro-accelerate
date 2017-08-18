#ifndef DEDISPERSE_KERNEL_H_
#define DEDISPERSE_KERNEL_H_

#define ARRAYSIZE SDIVINT * SDIVINDM

#include "float.h"

// Stores temporary shift values
__device__ __constant__ float dm_shifts[8192];
__device__ __constant__ int i_nsamp, i_nchans, i_t_processed_s;
__device__  __shared__ uchar4 f_line[UNROLLS][ARRAYSIZE + 1];

//{{{ shared_dedisperse_loop

__global__ void shared_dedisperse_kernel(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	unsigned char temp_f;

	int i, j, c, local_one, local_two, unroll, stage;

	int local_kernel_one[SNUMREG];
	int local_kernel_two[SNUMREG];
	int local_kernel_three[SNUMREG];
	int local_kernel_four[SNUMREG];

	int shift[UNROLLS];

	for (i = 0; i < SNUMREG; i++)
	{
		local_kernel_one[i] = 0;
		local_kernel_two[i] = 0;
		local_kernel_three[i] = 0;
		local_kernel_four[i] = 0;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT ) );

	float findex = ( threadIdx.x * 4) + 3;

	int nsamp_counter = ( idx + ( blockIdx.x * ( 4 * SNUMREG * SDIVINT ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS)
	{
		__syncthreads();

		for (j = 0; j < UNROLLS; j++)
		{
			temp_f = (unsigned char)( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) )  + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j][idx + 3].x = temp_f;
			f_line[j][idx + 2].y = temp_f;
			f_line[j][idx + 1].z = temp_f;
			f_line[j][idx    ].w = temp_f;

			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG; i++)
		{
			local_one = 0;
			local_two = 0;
			unroll = ( i * 4 * SDIVINT );
			for (j = 0; j < UNROLLS; j++)
			{
				stage = *(int*) &f_line[j][( shift[j] + unroll )];
				local_one += (stage & 0x00FF00FF);
				local_two += (stage & 0xFF00FF00)>>8;
			}
			local_kernel_one[i]   += (local_one & 0x0000FFFF);
			local_kernel_two[i]   += (local_one & 0xFFFF0000) >>16;
			local_kernel_three[i] += (local_two & 0x0000FFFF);
			local_kernel_four[i]  += (local_two & 0xFFFF0000) >> 16;
		}
	}

	// Write the accumulators to the output array. 
	local_one = ( ( ( ( blockIdx.y * SDIVINDM ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 4 * SNUMREG * SDIVINT ) ) + 4 * threadIdx.x;

	#pragma unroll
	for (i = 0; i < SNUMREG; i++)
	{
		*( (float2*) ( d_output + local_one + ( i * 4 * SDIVINT ) ) ) = make_float2((float)local_kernel_one[i] / i_nchans/bin,
											    (float)local_kernel_two[i] / i_nchans/bin);
		*( (float2*) ( d_output + local_one + 2 + ( i * 4 * SDIVINT ) ) ) = make_float2((float)local_kernel_three[i] / i_nchans/bin,
												(float)local_kernel_three[i] / i_nchans/bin);
	}
}


__global__ void shared_dedisperse_kernel_16(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	int i, c;
	int shift;

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG];
	float local_kernel_two[SNUMREG];

	for (i = 0; i < SNUMREG; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG * SDIVINT ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c ++)
	{

		__syncthreads();

			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c] * shift_two) ) ) + ( nsamp_counter )) );

			f_line[0][idx].x = temp_f;
			if (idx > 0)
			{
				f_line[0][idx - 1].y = temp_f;
			}
			shift = __float2int_rz(shift_one * dm_shifts[c] + findex);

		nsamp_counter = ( nsamp_counter + i_nsamp );

		__syncthreads();

		for (i = 0; i < SNUMREG; i++)
		{
			unroll = ( i * 2 * SDIVINT );
			local = *(int*) &f_line[0][( shift + unroll )];
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG * SDIVINT ) ) + 2 * threadIdx.x;

	#pragma unroll
	for (i = 0; i < SNUMREG; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
	}
}

__global__ void cache_dedisperse_kernel(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	int   shift;	
	float local_kernel;

	int t  = blockIdx.x * SDIVINT  + threadIdx.x;
	
	// Initialise the time accumulators
	local_kernel = 0.0f;

	float shift_temp = mstartdm + ((blockIdx.y * SDIVINDM + threadIdx.y) * mdmstep);
	
	// Loop over the frequency channels.
        for(int c = 0; c < i_nchans; c++) {


		// Calculate the initial shift for this given frequency
		// channel (c) at the current despersion measure (dm) 
		// ** dm is constant for this thread!!**
		shift = (c * (i_nsamp) + t) + __float2int_rz (dm_shifts[c] * shift_temp);
		
		local_kernel += (float)__ldg(&d_input[shift]);
	}

	// Write the accumulators to the output array. 
	shift = ( ( ( blockIdx.y * SDIVINDM ) + threadIdx.y ) * ( i_t_processed_s ) ) + t;

	d_output[shift] = (local_kernel / i_nchans / bin);

}
#endif



