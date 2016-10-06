#ifndef DEDISPERSE_KERNEL_H_
#define DEDISPERSE_KERNEL_H_

#define ARRAYSIZE SDIVINT * SDIVINDM

// Stores temporary shift values
__device__ __constant__ float dm_shifts[15500];
__device__ __constant__ int   i_nsamp, i_nchans, i_t_processed_s;
__device__ __shared__ float2 fa_line[ARRAYSIZE+1];
__device__ __shared__ float2 fb_line[ARRAYSIZE+1];
__device__ __shared__ float2 fc_line[ARRAYSIZE+1];
__device__ __shared__ float2 fd_line[ARRAYSIZE+1];

//{{{ shared_dedisperse_loop

__global__ void shared_dedisperse_kernel(unsigned char *d_input, float *d_output, cudaTextureObject_t tex, float mstartdm, float mdmstep)
{
//	int   i, c, shift;
	int   i, c, shifta, shiftb, shiftc, shiftd;

	float2 local_kernel_t[SNUMREG];

	#pragma unroll
	for(i = 0; i < SNUMREG; i++) local_kernel_t[i] = make_float2(0.0f, 0.0f);

	int idx = (threadIdx.x + (threadIdx.y * SDIVINT));
	int nsamp_counter = (idx + (blockIdx.x*(2*SNUMREG*SDIVINT)));

	float shift_two = (mstartdm + (__int2float_rz(blockIdx.y)*SFDIVINDM*mdmstep));
	float shift_one = (__int2float_rz(threadIdx.y)*mdmstep);

	for(c = 0; c < i_nchans; c+=4) {

		__syncthreads();
		
		float temp_f1 = (float)__ldg(d_input + nsamp_counter + __float2int_rz(dm_shifts[c]*shift_two));
		float temp_f2 = (float)__ldg(d_input + nsamp_counter + i_nsamp + __float2int_rz(dm_shifts[c+1]*shift_two));
		float temp_f3 = (float)__ldg(d_input + nsamp_counter + 2*i_nsamp + __float2int_rz(dm_shifts[c+2]*shift_two));
		float temp_f4 = (float)__ldg(d_input + nsamp_counter + 3*i_nsamp + __float2int_rz(dm_shifts[c+3]*shift_two));

		fa_line[idx].x = temp_f1;
                fb_line[idx].x = temp_f2; 
                fc_line[idx].x = temp_f3; 
                fd_line[idx].x = temp_f4; 
		if(idx > 0) {
			fa_line[idx-1].y = temp_f1;
			fb_line[idx-1].y = temp_f2;
			fc_line[idx-1].y = temp_f3;
			fd_line[idx-1].y = temp_f4;
		}
	
		shifta = __float2int_rz(shift_one*dm_shifts[c]) + threadIdx.x*2;
		shiftb = __float2int_rz(shift_one*dm_shifts[c+1]) + threadIdx.x*2;
		shiftc = __float2int_rz(shift_one*dm_shifts[c+2]) + threadIdx.x*2;
		shiftd = __float2int_rz(shift_one*dm_shifts[c+3]) + threadIdx.x*2;
		
		nsamp_counter += 4*i_nsamp;

		__syncthreads();

		#pragma unroll
		for(i = 0; i < SNUMREG; i++) {
			float2 local_fa = fa_line[shifta + (i*2*SDIVINT)];
			float2 local_fb = fb_line[shiftb + (i*2*SDIVINT)];
			float2 local_fc = fc_line[shiftc + (i*2*SDIVINT)];
			float2 local_fd = fd_line[shiftd + (i*2*SDIVINT)];
			local_kernel_t[i].x += (local_fa.x+local_fb.x+local_fc.x+local_fd.x);
			local_kernel_t[i].y += (local_fa.y+local_fb.y+local_fc.y+local_fd.y);

		}
	}

	// Write the accumulators to the output array. 
	shifta = ((blockIdx.y*SDIVINDM) + threadIdx.y)*(i_t_processed_s) + blockIdx.x*2*SNUMREG*SDIVINT + 2*threadIdx.x;

	#pragma unroll
	for(i = 0; i < SNUMREG; i++) {
		d_output[shifta + (i*2*SDIVINT)] = local_kernel_t[i].x/i_nchans;
		d_output[shifta + (i*2*SDIVINT) + 1] = local_kernel_t[i].y/i_nchans;
	}
	if(blockIdx.x==0 && threadIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) printf("\n%d %d %d", i_t_processed_s, i_nsamp, i_nchans);
}
#endif

