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
	int i, j, c;
	int shift[UNROLLS];

	unsigned char temp_f;
	long int local;
	int unroll;

	float findex = ( threadIdx.x * 4 );
	
	int local_kernel_one[SNUMREG];
	int local_kernel_two[SNUMREG];
	int local_kernel_three[SNUMREG];
	int local_kernel_four[SNUMREG];

	for (i = 0; i < SNUMREG; i++)
	{
		local_kernel_one[i] = 0;
		local_kernel_two[i] = 0;
		local_kernel_three[i] = 0;
		local_kernel_four[i] = 0;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 4 * SNUMREG * SDIVINT ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS)
	{
		__syncthreads();

		for (j = 0; j < UNROLLS; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j][idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j][idx - 1].y = temp_f;
			} else if (idx > 1) {
				f_line[j][idx - 2].z = temp_f;
			} else if (idx > 2) {
				f_line[j][idx - 3].w = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG; i++)
		{
			local = 0;
			unroll = ( i * 4 * SDIVINT );
			//local = *(int*) &f_line[j][( shift[j] + unroll )];

			for (j = 0; j < UNROLLS; j++)
			{
			//	local += *(long int*) &f_line[j][( shift[j] + unroll )];
			//}
			//local_kernel_one[i]   += ( (ushort4*) ( &local ) )->x;
			//local_kernel_two[i]   += ( (ushort4*) ( &local ) )->y;
			//local_kernel_three[i] += ( (ushort4*) ( &local ) )->z;
			//local_kernel_four[i]  += ( (ushort4*) ( &local ) )->w;


			local = *(int*) &f_line[j][( shift[j] + unroll )];
			local_kernel_one[i]   += ( (uchar4*) ( &local ) )->x;
			local_kernel_two[i]   += ( (uchar4*) ( &local ) )->y;
			local_kernel_three[i] += ( (uchar4*) ( &local ) )->z;
			local_kernel_four[i]  += ( (uchar4*) ( &local ) )->w;
			}
		}
	}

	// Write the accumulators to the output array. 
	int loc = ( ( ( ( blockIdx.y * SDIVINDM ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 4 * SNUMREG * SDIVINT ) ) + 4 * threadIdx.x;

	#pragma unroll
	for (i = 0; i < SNUMREG; i++)
	{
		*( (float4*) ( d_output + loc + ( i * 4 * SDIVINT ) ) ) = make_float4(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin, local_kernel_three[i] / i_nchans / bin, local_kernel_four[i] / i_nchans / bin);
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



