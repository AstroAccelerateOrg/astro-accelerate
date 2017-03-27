#ifndef DEDISPERSE_KERNEL_H_
#define DEDISPERSE_KERNEL_H_

#define ARRAYSIZE SDIVINT * SDIVINDM

#include "float.h"
#include "headers/kernel_params.h"

// Stores temporary shift values
__device__ __constant__ float dm_shifts[8192];
__device__ __constant__ int i_nsamp, i_nchans, i_t_processed_s;
//__device__  __shared__ ushort2 f_line[UNROLLS][ARRAYSIZE + 1];

//{{{ shared_dedisperse_loop

__global__ void shared_dedisperse_kernel_range_0(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS_0];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG_0];
	float local_kernel_two[SNUMREG_0];

	for (i = 0; i < SNUMREG_0; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT_0 ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG_0 * SDIVINT_0 ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM_0 * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS_0)
	{

		__syncthreads();

		for (j = 0; j < UNROLLS_0; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS_0 + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS_0 + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS_0 * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG_0; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT_0 );
			for (j = 0; j < UNROLLS_0; j++)
			{
				local += *(int*) &f_line[j*UNROLLS_0 + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM_0 ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG_0 * SDIVINT_0 ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG_0; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT_0 ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT_0)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT_0) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


__global__ void shared_dedisperse_kernel_range_1(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS_1];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG_1];
	float local_kernel_two[SNUMREG_1];

	for (i = 0; i < SNUMREG_1; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT_1 ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG_1 * SDIVINT_1 ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM_1 * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS_1)
	{

		__syncthreads();

		for (j = 0; j < UNROLLS_1; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS_1 + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS_1 + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS_1 * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG_1; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT_1 );
			for (j = 0; j < UNROLLS_1; j++)
			{
				local += *(int*) &f_line[j*UNROLLS_1 + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM_1 ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG_1 * SDIVINT_1 ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG_1; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT_1 ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT_1)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT_1) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


__global__ void shared_dedisperse_kernel_range_2(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS_2];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG_2];
	float local_kernel_two[SNUMREG_2];

	for (i = 0; i < SNUMREG_2; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT_2 ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG_2 * SDIVINT_2 ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM_2 * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS_2)
	{

		__syncthreads();

		for (j = 0; j < UNROLLS_2; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS_2 + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS_2 + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS_2 * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG_2; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT_2 );
			for (j = 0; j < UNROLLS_2; j++)
			{
				local += *(int*) &f_line[j*UNROLLS_2 + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM_2 ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG_2 * SDIVINT_2 ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG_2; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT_2 ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT_2)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT_2) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


__global__ void shared_dedisperse_kernel_range_3(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS_3];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG_3];
	float local_kernel_two[SNUMREG_3];

	for (i = 0; i < SNUMREG_3; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT_3 ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG_3 * SDIVINT_3 ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM_3 * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS_3)
	{

		__syncthreads();

		for (j = 0; j < UNROLLS_3; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS_3 + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS_3 + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS_3 * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG_3; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT_3 );
			for (j = 0; j < UNROLLS_3; j++)
			{
				local += *(int*) &f_line[j*UNROLLS_3 + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM_3 ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG_3 * SDIVINT_3 ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG_3; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT_3 ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT_3)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT_3) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


__global__ void shared_dedisperse_kernel_range_4(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS_4];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG_4];
	float local_kernel_two[SNUMREG_4];

	for (i = 0; i < SNUMREG_4; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT_4 ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG_4 * SDIVINT_4 ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM_4 * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS_4)
	{

		__syncthreads();

		for (j = 0; j < UNROLLS_4; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS_4 + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS_4 + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS_4 * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG_4; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT_4 );
			for (j = 0; j < UNROLLS_4; j++)
			{
				local += *(int*) &f_line[j*UNROLLS_4 + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM_4 ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG_4 * SDIVINT_4 ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG_4; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT_4 ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT_4)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT_4) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


__global__ void shared_dedisperse_kernel_range_5(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS_5];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG_5];
	float local_kernel_two[SNUMREG_5];

	for (i = 0; i < SNUMREG_5; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT_5 ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG_5 * SDIVINT_5 ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM_5 * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS_5)
	{

		__syncthreads();

		for (j = 0; j < UNROLLS_5; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS_5 + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS_5 + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS_5 * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG_5; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT_5 );
			for (j = 0; j < UNROLLS_5; j++)
			{
				local += *(int*) &f_line[j*UNROLLS_5 + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM_5 ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG_5 * SDIVINT_5 ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG_5; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT_5 ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT_5)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT_5) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


__global__ void shared_dedisperse_kernel_range_6(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS_6];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG_6];
	float local_kernel_two[SNUMREG_6];

	for (i = 0; i < SNUMREG_6; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT_6 ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG_6 * SDIVINT_6 ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM_6 * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS_6)
	{

		__syncthreads();

		for (j = 0; j < UNROLLS_6; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS_6 + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS_6 + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS_6 * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG_6; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT_6 );
			for (j = 0; j < UNROLLS_6; j++)
			{
				local += *(int*) &f_line[j*UNROLLS_6 + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM_6 ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG_6 * SDIVINT_6 ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG_6; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT_6 ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT_6)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT_6) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


__global__ void shared_dedisperse_kernel_range_7(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS_7];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG_7];
	float local_kernel_two[SNUMREG_7];

	for (i = 0; i < SNUMREG_7; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT_7 ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG_7 * SDIVINT_7 ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM_7 * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS_7)
	{

		__syncthreads();

		for (j = 0; j < UNROLLS_7; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS_7 + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS_7 + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS_7 * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG_7; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT_7 );
			for (j = 0; j < UNROLLS_7; j++)
			{
				local += *(int*) &f_line[j*UNROLLS_7 + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM_7 ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG_7 * SDIVINT_7 ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG_7; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT_7 ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT_7)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT_7) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


#endif

