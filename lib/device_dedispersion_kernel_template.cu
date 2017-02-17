__global__ void shared_dedisperse_kernel_range{RANGE}(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS{RANGE}];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG{RANGE}];
	float local_kernel_two[SNUMREG{RANGE}];

	for (i = 0; i < SNUMREG{RANGE}; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT{RANGE} ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG{RANGE} * SDIVINT{RANGE} ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM{RANGE} * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS{RANGE})
	{

		__syncthreads();

		for (j = 0; j < UNROLLS{RANGE}; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS{RANGE} + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS{RANGE} + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS{RANGE} * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG{RANGE}; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT{RANGE} );
			for (j = 0; j < UNROLLS{RANGE}; j++)
			{
				local += *(int*) &f_line[j*UNROLLS{RANGE} + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM{RANGE} ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG{RANGE} * SDIVINT{RANGE} ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG{RANGE}; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT{RANGE} ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT{RANGE})    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT{RANGE}) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


