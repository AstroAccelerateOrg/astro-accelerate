#include "aa_device_dedispersion_kernel.hpp"

#include "float.h"
#include <stdio.h>
#include "aa_timelog.hpp"

namespace astroaccelerate {

  //{{{ shared_dedisperse_loop

  __device__ __shared__ ushort2 f_line[UNROLLS][ARRAYSIZE + 2];
  __device__ __shared__ uchar4 test[UNROLLS][ARRAYSIZE + 4];
  __device__ __constant__ int i_nsamp, i_nchans, i_t_processed_s;
  __device__ __constant__ float dm_shifts[8192];

  __global__ void shared_dedisperse_kernel_4bit_4096p(unsigned short *d_input, float *d_output, float *d_dm_shifts, float mstartdm, float mdmstep, int chan_shift){

     unsigned char temp_f;

        int i, j, c, shift_c, local, unroll, stage;

        int shift[UNROLLS];
        int local_kernel_one[SNUMREG];
        int local_kernel_two[SNUMREG];
	unsigned int sum_results[4][SNUMREG];

        float findex = ((threadIdx.x*4) + 3); // + 1 

        int idx           = (threadIdx.x + (threadIdx.y*SDIVINT));
        int nsamp_counter = (idx + (blockIdx.x*(4*SNUMREG*SDIVINT)));

        float shift_two = (mstartdm + (__int2float_rz(blockIdx.y)*SFDIVINDM*mdmstep));
        float shift_one = (__int2float_rz(threadIdx.y)*mdmstep);

	for (i = 0; i < SNUMREG; i++){
		sum_results[0][i] = 0;
		sum_results[1][i] = 0;
		sum_results[2][i] = 0;
		sum_results[3][i] = 0;
	}

	for (shift_c = 0; shift_c < chan_shift; shift_c++){
	        for (i = 0; i < SNUMREG; i++){
	                local_kernel_one[i] = 0;
	                local_kernel_two[i] = 0;
	        }
        for (c = shift_c*4096; c < (shift_c + 1)*4096 ; c += UNROLLS)
        {
		if (c < i_nchans){
                __syncthreads();

                for (j = 0; j < UNROLLS; j++)
                {
                        temp_f = (unsigned char)(__ldg((d_input + (__float2int_rz(d_dm_shifts[c + j]*shift_two))) + (nsamp_counter + (j*i_nsamp))));

                        test[j][idx + 3].x = temp_f;
                        test[j][idx + 2].y = temp_f;
                        test[j][idx + 1].z = temp_f;
                        test[j][idx    ].w = temp_f;

                        shift[j] = __float2int_rz(shift_one*d_dm_shifts[c + j] + findex);
                }

                nsamp_counter = (nsamp_counter + (UNROLLS*i_nsamp));

                __syncthreads();

                for (i = 0; i < SNUMREG; i++)
                {
                        local = 0;
                        unroll = (i*4*SDIVINT);
                        for (j = 0; j < UNROLLS; j++)
                        {
                                stage = *(int*) &test[j][(shift[j] + unroll)];
                                local += stage;
                        }
                        local_kernel_one[i] += (local & 0x00FF00FF);
                        local_kernel_two[i] += ((local & 0xFF00FF00) >> 8);
                }
		}
        }

	for (i = 0; i < SNUMREG; i++){
		sum_results[0][i] += (local_kernel_one[i] &0x0000FFFF);
		sum_results[1][i] += (local_kernel_two[i] &0x0000FFFF);
		sum_results[2][i] += (local_kernel_one[i] &0xFFFF0000) >> 16;
		sum_results[3][i] += (local_kernel_two[i] &0xFFFF0000) >> 16;
	}
	} // shift in chan


        // Write the accumulators to the output array.
        local = ((((blockIdx.y*SDIVINDM ) + threadIdx.y)*(i_t_processed_s)) + (blockIdx.x*4*SNUMREG*SDIVINT)) + 4*threadIdx.x;

        #pragma unroll
        for (i = 0; i < SNUMREG; i++){
                *((float4*)(d_output + local + (i*4*SDIVINT ))) = make_float4((float)(sum_results[0][i])/i_nchans,
                                                                              (float)(sum_results[1][i])/i_nchans,
                                                                              (float)(sum_results[2][i])/i_nchans,
                                                                              (float)(sum_results[3][i])/i_nchans);
        }

}

  __global__ void shared_dedisperse_kernel_4bit(unsigned short *d_input, float *d_output, float mstartdm, float mdmstep){

     unsigned char temp_f;

        int i, j, c, local, unroll, stage;

        int shift[UNROLLS];
        int local_kernel_one[SNUMREG];
        int local_kernel_two[SNUMREG];

        float findex = (( threadIdx.x * 4 ) + 3 ); // + 1 

        for (i = 0; i < SNUMREG; i++)
        {
                local_kernel_one[i] = 0;
                local_kernel_two[i] = 0;
        }

        int idx           = (threadIdx.x + (threadIdx.y*SDIVINT));
        int nsamp_counter = (idx + (blockIdx.x*(4*SNUMREG*SDIVINT)));

        float shift_two = (mstartdm + (__int2float_rz(blockIdx.y)*SFDIVINDM*mdmstep));
        float shift_one = (__int2float_rz(threadIdx.y)*mdmstep);

        for (c = 0; c < i_nchans; c += UNROLLS)
        {
                __syncthreads();

                for (j = 0; j < UNROLLS; j++)
                {
                        temp_f = (unsigned char)(__ldg((d_input + (__float2int_rz(dm_shifts[c + j]*shift_two))) + (nsamp_counter + (j*i_nsamp))));

                        test[j][idx + 3].x = temp_f;
                        test[j][idx + 2].y = temp_f;
                        test[j][idx + 1].z = temp_f;
                        test[j][idx    ].w = temp_f;

                        shift[j] = __float2int_rz(shift_one*dm_shifts[c + j] + findex);
                }

                nsamp_counter = (nsamp_counter + (UNROLLS*i_nsamp));

                __syncthreads();

                for (i = 0; i < SNUMREG; i++)
                {
                        local = 0;
                        unroll = (i*4*SDIVINT);
                        for (j = 0; j < UNROLLS; j++)
                        {
                                stage = *(int*) &test[j][(shift[j] + unroll)];
                                local += stage;
                        }
                        local_kernel_one[i] += (local & 0x00FF00FF);
                        local_kernel_two[i] += ((local & 0xFF00FF00) >> 8);
                }
        }

        // Write the accumulators to the output array.
        local = ((((blockIdx.y*SDIVINDM) + threadIdx.y)*(i_t_processed_s)) + (blockIdx.x*4*SNUMREG*SDIVINT)) + 4 * threadIdx.x;

        #pragma unroll
        for (i = 0; i < SNUMREG; i++)
        {
                *( (float4*) ( d_output + local + ( i * 4 * SDIVINT ) ) ) = make_float4((float)(local_kernel_one[i] &0x0000FFFF) / i_nchans,
                                                                                        (float)( (local_kernel_two[i] &0x0000FFFF)) / i_nchans,
                                                                                        (float)( (local_kernel_one[i] &0xFFFF0000) >> 16) / i_nchans,
                                                                                        (float)( (local_kernel_two[i] &0xFFFF0000) >> 16) / i_nchans);
        }
}

  __global__ void shared_dedisperse_kernel(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep) {
    ushort temp_f;

    int i, j, c, local, unroll, stage;

    int shift[UNROLLS];
    int local_kernel_one[SNUMREG];
    int local_kernel_two[SNUMREG];

    float findex = (( threadIdx.x * 2 ) + 1 );

    for (i = 0; i < SNUMREG; i++)
      {
	local_kernel_one[i] = 0;
	local_kernel_two[i] = 0;
      }

    int idx 	  = ( threadIdx.x + ( threadIdx.y * SDIVINT ) );
    int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG * SDIVINT ) ) );

    float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM * mdmstep ) );
    float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );

    for (c = 0; c < i_nchans; c += UNROLLS)
      {

	__syncthreads();

	for (j = 0; j < UNROLLS; j++)
	  {
	    temp_f = ( __ldg(( d_input + ( __float2int_rz(dm_shifts[c + j] * shift_two) ) )  + ( nsamp_counter + ( j * i_nsamp ) )) );

	    f_line[j][idx + 1].x = temp_f;
	    f_line[j][idx    ].y = temp_f;

	    shift[j] = __float2int_rz(shift_one * dm_shifts[c + j] + findex);
	  }

	nsamp_counter = ( nsamp_counter + ( UNROLLS * i_nsamp ) );

	__syncthreads();

	for (i = 0; i < SNUMREG; i++)
	  {
	    local = 0;
	    unroll = ( i * 2 * SDIVINT );
	    for (j = 0; j < UNROLLS; j++)
	      {
		stage = *(int*) &f_line[j][( shift[j] + unroll )];
		local += stage;
	      }
	    local_kernel_one[i] += (local & 0x0000FFFF);
	    local_kernel_two[i] += (local & 0xFFFF0000) >> 16;
	  }
      }

    // Write the accumulators to the output array. 
    size_t big_local;
    big_local = ( ( (size_t)( ( blockIdx.y * SDIVINDM ) + threadIdx.y ) * ( (size_t)i_t_processed_s ) ) + (size_t)( blockIdx.x * 2 * SNUMREG * SDIVINT ) ) + (size_t)(2 * threadIdx.x);

#pragma unroll
    for (i = 0; i < SNUMREG; i++)
      {
	*( (float2*) ( d_output + big_local + (size_t)( i*2*SDIVINT ) ) ) = make_float2((float)local_kernel_one[i] / i_nchans/bin,(float)local_kernel_two[i] / i_nchans/bin);
      }
  }

	__global__ void shared_dedisperse_kernel_nchan8192p(int bin, unsigned short *d_input, float *d_output, float *d_dm_shifts, float mstartdm, float mdmstep) {
		ushort temp_f;
		
		int i, j, c, local, unroll, stage;
		
		int shift[UNROLLS];
		int local_kernel_one[SNUMREG];
		int local_kernel_two[SNUMREG];
		
		float findex = (( threadIdx.x * 2 ) + 1 );
		
		for (i = 0; i < SNUMREG; i++) {
			local_kernel_one[i] = 0;
			local_kernel_two[i] = 0;
		}
		
		int idx 	  = ( threadIdx.x + ( threadIdx.y * SDIVINT ) );
		int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG * SDIVINT ) ) );
		
		float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM * mdmstep ) );
		float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );
		
		for (c = 0; c < i_nchans; c += UNROLLS) {
			
			__syncthreads();
			for (j = 0; j < UNROLLS; j++) {
				temp_f = ( __ldg(( d_input + ( __float2int_rz(d_dm_shifts[c + j] * shift_two) ) )  + ( nsamp_counter + ( j * i_nsamp ) )) );
				
				f_line[j][idx + 1].x = temp_f;
				f_line[j][idx    ].y = temp_f;
				
				shift[j] = __float2int_rz(shift_one * d_dm_shifts[c + j] + findex);
			}
			
			nsamp_counter = ( nsamp_counter + ( UNROLLS * i_nsamp ) );
			
			__syncthreads();
		
			for (i = 0; i < SNUMREG; i++) {
				local = 0;
				unroll = ( i * 2 * SDIVINT );
				for (j = 0; j < UNROLLS; j++) {
					stage = *(int*) &f_line[j][( shift[j] + unroll )];
					local += stage;
				}
				local_kernel_one[i] += (local & 0x0000FFFF);
				local_kernel_two[i] += (local & 0xFFFF0000) >> 16;
			}
		}
		
		// Write the accumulators to the output array. 
		size_t big_local;
		big_local = (((size_t)((blockIdx.y*SDIVINDM ) + threadIdx.y)*((size_t)i_t_processed_s)) + (size_t)(blockIdx.x*2*SNUMREG*SDIVINT)) + (size_t)(2*threadIdx.x);
		
		#pragma unroll
		for (i = 0; i < SNUMREG; i++) {
			*( (float2*) ( d_output + big_local + (size_t)(i*2*SDIVINT) ) ) = make_float2((float)local_kernel_one[i]/i_nchans/bin, (float)local_kernel_two[i]/i_nchans/bin);
		}
	}

  __global__ void shared_dedisperse_kernel_16(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep) {
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
	size_t big_local;
    big_local = ( (size_t)( ( ( blockIdx.y * SDIVINDM ) + threadIdx.y )*( i_t_processed_s ) ) + (size_t)( blockIdx.x * 2 * SNUMREG * SDIVINT ) ) + (size_t)(2*threadIdx.x);

#pragma unroll
    for (i = 0; i < SNUMREG; i++)
      {
	*( (float2*) ( d_output + big_local + (size_t)(i*2*SDIVINT ) ) ) = make_float2(local_kernel_one[i] / i_nchans / bin, local_kernel_two[i] / i_nchans / bin);
      }
  }

  
	__global__ void shared_dedisperse_kernel_16_nchan8192p(int bin, unsigned short *d_input, float *d_output, float *d_dm_shifts, float mstartdm, float mdmstep) {
		int i, c;
		int shift;
		
		ushort temp_f;
		int local, unroll;
		
		float findex = ( threadIdx.x * 2 );
		float local_kernel_one[SNUMREG];
		float local_kernel_two[SNUMREG];
		
		for (i = 0; i < SNUMREG; i++) {
			local_kernel_one[i] = 0.0f;
			local_kernel_two[i] = 0.0f;
		}
		
		int idx = ( threadIdx.x + ( threadIdx.y * SDIVINT ) );
		int nsamp_counter = ( idx + ( blockIdx.x * ( 2 * SNUMREG * SDIVINT ) ) );
		
		float shift_two = ( mstartdm + ( __int2float_rz(blockIdx.y) * SFDIVINDM * mdmstep ) );
		float shift_one = ( __int2float_rz(threadIdx.y) * mdmstep );
		
		for (c = 0; c < i_nchans; c ++) {
			
			__syncthreads();
			
			temp_f = ( __ldg(( d_input + ( __float2int_rz(d_dm_shifts[c] * shift_two) ) ) + ( nsamp_counter )) );
			
			f_line[0][idx].x = temp_f;
			if (idx > 0) {
				f_line[0][idx - 1].y = temp_f;
			}
			shift = __float2int_rz(shift_one * d_dm_shifts[c] + findex);
			
			nsamp_counter = ( nsamp_counter + i_nsamp );
			
			__syncthreads();
			
			for (i = 0; i < SNUMREG; i++) {
				unroll = ( i * 2 * SDIVINT );
				local = *(int*) &f_line[0][( shift + unroll )];
				local_kernel_one[i] += ( (ushort2*) ( &local ) )->x;
				local_kernel_two[i] += ( (ushort2*) ( &local ) )->y;
			}
		}
		
		// Write the accumulators to the output array. 
		size_t big_local;
		big_local = ( (size_t)( ( ( blockIdx.y * SDIVINDM ) + threadIdx.y ) * ( i_t_processed_s ) ) + (size_t)( blockIdx.x * 2 * SNUMREG * SDIVINT ) ) + (size_t)(2*threadIdx.x);
		
		#pragma unroll
		for (i = 0; i < SNUMREG; i++) {
			*( (float2*) ( d_output + big_local + (size_t)( i * 2 * SDIVINT ) ) ) = make_float2(local_kernel_one[i]/i_nchans/bin, local_kernel_two[i]/i_nchans/bin);
		}
	}

  
  __global__ void cache_dedisperse_kernel(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep) {
    size_t   shift;	
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

      shift = ((size_t)(c*(i_nsamp)) + (size_t)t) + (size_t)(__float2int_rz (dm_shifts[c] * shift_temp));
		
      local_kernel += (float)__ldg(&d_input[shift]);
    }

    // Write the accumulators to the output array. 
    shift = ( ( (size_t)( blockIdx.y * SDIVINDM ) + (size_t)threadIdx.y ) * ( (size_t)i_t_processed_s ) ) + (size_t)t;

    d_output[shift] = (local_kernel / i_nchans / bin);

  }
  
	__global__ void cache_dedisperse_kernel_nchan8192p(int bin, unsigned short *d_input, float *d_output, float *d_dm_shifts, float mstartdm, float mdmstep) {
		size_t   shift;	
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
			shift = ((size_t)(c*(i_nsamp)) + (size_t)t) + (size_t)(__float2int_rz (d_dm_shifts[c]*shift_temp));
			
			local_kernel += (float)__ldg(&d_input[shift]);
		}

		// Write the accumulators to the output array. 
		shift = ( (size_t)( ( blockIdx.y * SDIVINDM ) + threadIdx.y ) * ( (size_t)i_t_processed_s ) ) + (size_t)t;

		d_output[shift] = (local_kernel / i_nchans / bin);
  }
 

	//-------------------------------- wrapper functions
  
  /** \brief Kernel wrapper function to set device constants for dedispersion_kernel kernel function. */
  void set_device_constants_dedispersion_kernel(const int &nchans, const int &length, const int &t_processed, const float *const dmshifts) {
    cudaMemcpyToSymbol(dm_shifts, dmshifts, nchans * sizeof(float));
    cudaMemcpyToSymbol(i_nchans, &nchans, sizeof(int));
    cudaMemcpyToSymbol(i_nsamp, &length, sizeof(int));
    cudaMemcpyToSymbol(i_t_processed_s, &t_processed, sizeof(int));
    //checkCudaErrors(cudaGetLastError());
  }

  /** \brief Kernel wrapper function to set device constants for kernel dedispersion_kernel function. */
  void set_device_constants_dedispersion_kernel(const long int &length, const int &t_processed) {
    cudaMemcpyToSymbol(i_nsamp, &length, sizeof(int));
    cudaMemcpyToSymbol(i_t_processed_s, &t_processed, sizeof(int));
  }

  void set_device_constants_dedispersion_kernel(const int &nchans, const long int &length, const int &t_processed) {
    cudaMemcpyToSymbol(i_nchans, &nchans, sizeof(int));
    cudaMemcpyToSymbol(i_nsamp, &length, sizeof(int));
    cudaMemcpyToSymbol(i_t_processed_s, &t_processed, sizeof(int));
  }

  /** \brief Kernel wrapper function for dedisperse_kernel  kernel function. */
  void call_kernel_shared_dedisperse_kernel(const dim3 &block_size, const dim3 &grid_size,
					    const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep) {
    cudaFuncSetCacheConfig(shared_dedisperse_kernel, cudaFuncCachePreferShared);
    shared_dedisperse_kernel<<<block_size, grid_size>>>(bin, d_input, d_output, mstartdm, mdmstep);
  }

	/** \brief Kernel wrapper function for dedispersion GPU kernel which works with 4-bit input data */
	void call_kernel_shared_dedisperse_kernel_4bit(const dim3 &block_size, const dim3 &grid_size,
						  unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep){
		cudaFuncSetCacheConfig(shared_dedisperse_kernel_4bit, cudaFuncCachePreferShared);
		shared_dedisperse_kernel_4bit<<<block_size, grid_size>>>(d_input, d_output, mstartdm, mdmstep);
	}

        /** \brief Kernel wrapper function for dedispersion GPU kernel which works with 4-bit input data */
        void call_kernel_shared_dedisperse_kernel_4bit_4096chan(const dim3 &block_size, const dim3 &grid_size,
                                                  unsigned short *const d_input, float *const d_output, float *const d_dm_shifts, const float &mstartdm, const float &mdmstep, const int nchans){
		unsigned int chan_counter = (int)(ceil(nchans/4096));
                cudaFuncSetCacheConfig(shared_dedisperse_kernel_4bit_4096p, cudaFuncCachePreferShared);
                shared_dedisperse_kernel_4bit_4096p<<<block_size, grid_size>>>(d_input, d_output, d_dm_shifts, mstartdm, mdmstep, chan_counter);
        }
  
	/** \brief Kernel wrapper function for dedispersion GPU kernel which works with number of channels greater than 8192. */
	void call_kernel_shared_dedisperse_kernel_nchan8192p(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, float *const d_dm_shifts, const float &mstartdm, const float &mdmstep) {
		cudaFuncSetCacheConfig(shared_dedisperse_kernel_nchan8192p, cudaFuncCachePreferShared);
		shared_dedisperse_kernel_nchan8192p<<<block_size, grid_size>>>(bin, d_input, d_output, d_dm_shifts, mstartdm, mdmstep);
	}
  

  /** \brief Kernel wrapper function for dedisperse_kernel_16 kernel function. */
  void call_kernel_shared_dedisperse_kernel_16(const dim3 &block_size, const dim3 &grid_size,
					       const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep) {
    cudaFuncSetCacheConfig(shared_dedisperse_kernel_16, cudaFuncCachePreferShared);
    shared_dedisperse_kernel_16<<<block_size, grid_size>>>(bin, d_input, d_output, mstartdm, mdmstep);
  }
  
	/** \brief Kernel wrapper function for dedispersion kernel which works with 16bit data and when number of channels is greater than 8192 kernel function. */
	void call_kernel_shared_dedisperse_kernel_16_nchan8192p(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, float *const d_dm_shifts, const float &mstartdm, const float &mdmstep) {
		cudaFuncSetCacheConfig(shared_dedisperse_kernel_16_nchan8192p, cudaFuncCachePreferShared);
		shared_dedisperse_kernel_16<<<block_size, grid_size>>>(bin, d_input, d_output, mstartdm, mdmstep);
	}

  /** \brief Kernel wrapper function for cache_dedisperse_kernel kernel function. */
  void call_kernel_cache_dedisperse_kernel(const dim3 &block_size, const dim3 &grid_size,
					   const int &bin, unsigned short *const d_input, float *const d_output, const float &mstartdm, const float &mdmstep) {
    cudaFuncSetCacheConfig(cache_dedisperse_kernel, cudaFuncCachePreferL1);
    cache_dedisperse_kernel<<<block_size, grid_size>>>(bin, d_input, d_output, mstartdm, mdmstep);
  }
  
	/** \brief Kernel wrapper function for cache_dedisperse_kernel kernel function. */
	void call_kernel_cache_dedisperse_kernel_nchan8192p(const dim3 &block_size, const dim3 &grid_size, const int &bin, unsigned short *const d_input, float *const d_output, float *const d_dm_shifts, const float &mstartdm, const float &mdmstep) {
		cudaFuncSetCacheConfig(cache_dedisperse_kernel_nchan8192p, cudaFuncCachePreferL1);
		cache_dedisperse_kernel_nchan8192p<<<block_size, grid_size>>>(bin, d_input, d_output, d_dm_shifts, mstartdm, mdmstep);
	}

} //namespace astroaccelerate
