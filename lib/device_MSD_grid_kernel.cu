// Added by Karel Adamek 
#ifndef MSD_GRID_KERNEL_H_
#define MSD_GRID_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

__global__ void MSD_GPU_grid(float const* __restrict__ d_input, float *d_output, int x_steps, int y_steps, int nColumns, int msd)
{
	extern __shared__ float Ms_Ss[];

	int warp_id, local_id, dim_y, pos;
	float x; // current element
	float M; // streaming mean
	float S; // streaming sum of squares (x_i-\bar{x})
	float j;
	float ftemp;

	local_id = threadIdx.x & ( WARP - 1 );
	warp_id = threadIdx.x >> 5;
	dim_y = blockDim.x >> 5;

	//                           y                              +                      x
	pos = ( blockIdx.y * dim_y + warp_id ) * y_steps * nColumns + blockIdx.x * WARP * x_steps + local_id;
	M = __ldg(&d_input[pos]);
	S = 0;
	j = 1.0f;
	for (int xf = 1; xf < x_steps; xf++)
	{
		pos = pos + WARP;
		x = __ldg(&d_input[pos]);
		j = j + 1.0f;
		M = M + x;
		ftemp = ( j * x - M );
		S = S + 1.0f / ( j * ( j - 1.0f ) ) * ftemp * ftemp;
	}

	pos = pos + nColumns - ( x_steps - 1 ) * WARP;
	for (int yf = 1; yf < y_steps; yf++)
	{
		for (int xf = 0; xf < x_steps; xf++)
		{
			x = __ldg(&d_input[pos]);
			j = j + 1.0f;
			M = M + x;
			ftemp = ( j * x - M );
			S = S + 1.0f / ( j * ( j - 1.0f ) ) * ftemp * ftemp;
			pos = pos + WARP;
		}
		pos = pos + nColumns - x_steps * WARP;
	}

	Ms_Ss[threadIdx.x] = M;
	Ms_Ss[blockDim.x + threadIdx.x] = S;

	__syncthreads();

	// now all threads had saved their work, reduction follows

	for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1)
	{
		if (threadIdx.x < i)
		{
			j = j * 2;
			ftemp = ( M - Ms_Ss[i + threadIdx.x] );
			S = S + Ms_Ss[blockDim.x + i + threadIdx.x] + ( 1.0f / j ) * ftemp * ftemp;
			M = M + Ms_Ss[i + threadIdx.x];

			Ms_Ss[threadIdx.x] = M;
			Ms_Ss[blockDim.x + threadIdx.x] = S;
		}
		__syncthreads();
	}

	// by now we should have only 32 partial results. Shuffle reduction follows
	for (int q = HALF_WARP; q > 0; q = q >> 1)
	{
		j = j * 2;
		ftemp = ( M - __shfl_down(M, q) );
		S = S + __shfl_down(S, q) + ( 1.0f / j ) * ftemp * ftemp;
		M = M + __shfl_down(M, q);
	}

	//----------------------------------------------
	//---- Writing data
	if (threadIdx.x == 0)
	{
		pos = blockIdx.y * gridDim.x + blockIdx.x;
		if (msd)
		{
			// produce mean and sd instead of T and S.
			d_output[2 * pos] = M / j;
			d_output[2 * pos + 1] = sqrt(S / j);
		}
		else
		{
			d_output[2 * pos] = M;
			d_output[2 * pos + 1] = S;
		}
	}
}

#endif

