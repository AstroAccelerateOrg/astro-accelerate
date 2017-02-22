__global__ void shared_dedisperse_kernel_range__RANGE__(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep)
{
	extern __shared__ ushort2 f_line[];
	int i, j, c;
	int shift[UNROLLS__RANGE__];

	ushort temp_f;
	int local, unroll;

	float findex = ( threadIdx.x * 2 );
	float local_kernel_one[SNUMREG__RANGE__];
	float local_kernel_two[SNUMREG__RANGE__];

	for (i = 0; i < SNUMREG__RANGE__; i++)
	{
		local_kernel_one[i] = 0.0f;
		local_kernel_two[i] = 0.0f;
	}

	int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT__RANGE__ ) );
	int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG__RANGE__ * SDIVINT__RANGE__ ) ) );

	float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM__RANGE__ * mdmstep ) );
	float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

	for (c = 0; c < i_nchans; c += UNROLLS__RANGE__)
	{

		__syncthreads();

		for (j = 0; j < UNROLLS__RANGE__; j++)
		{
			temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) ) + ( nsamp_counter + ( j * i_nsamp ) )) );

			f_line[j*UNROLLS__RANGE__ + idx].x = temp_f;
			if (idx > 0)
			{
				f_line[j*UNROLLS__RANGE__ + idx - 1].y = temp_f;
			}
			shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
		}

		nsamp_counter = ( nsamp_counter + ( UNROLLS__RANGE__ * i_nsamp ) );

		__syncthreads();

		for (i = 0; i < SNUMREG__RANGE__; i++)
		{
			local = 0;
			unroll = ( i * 2 * SDIVINT__RANGE__ );
			for (j = 0; j < UNROLLS__RANGE__; j++)
			{
				local += *(int*) &f_line[j*UNROLLS__RANGE__ + ( shift[j] + unroll )];
			}
			local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
			local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
		}
	}

	// Write the accumulators to the output array. 
	local = ( ( ( ( blockIdx.y * SDIVINDM__RANGE__ ) + threadIdx.y ) * ( i_t_processed_s ) ) + ( blockIdx.x * 2 * SNUMREG__RANGE__ * SDIVINT__RANGE__ ) ) + 2 * threadIdx.x;

#pragma unroll
	for (i = 0; i < SNUMREG__RANGE__; i++)
	{
		*( (float2*) ( d_output + local + ( i * 2 * SDIVINT__RANGE__ ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
//		d_output[local + (i*2*SDIVINT__RANGE__)    ] = (local_kernel_one[i])/i_nchans;
//		d_output[local + (i*2*SDIVINT__RANGE__) + 1] = (local_kernel_two[i])/i_nchans;
	}
}


