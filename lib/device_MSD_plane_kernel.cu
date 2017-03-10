// Added by Karel Adamek 

#ifndef MSD_PLANE_KERNEL_H_
#define MSD_PLANE_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

__global__ void MSD_GPU(float2 const* __restrict__ d_input, float *d_output)
{
	__shared__ float Ms[WARP * MSD_WARPS_PER_BLOCK];
	__shared__ float Ss[WARP * MSD_WARPS_PER_BLOCK];

	int warp_id, local_id, pos;
	float2 x;
	float M;
	float S;
	float j;
	float ftemp;

	local_id = threadIdx.x & ( WARP - 1 );
	warp_id = threadIdx.x >> 5;

	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos = blockIdx.x * MSD_WARPS_PER_BLOCK * WARP * MSD_ELEM_PER_THREAD + warp_id * WARP * MSD_ELEM_PER_THREAD + local_id;
	x = __ldg(&d_input[pos]);
	M = x.x;
	S = 0;
	j = 1.0f;

	j = j + 1.0f;
	M = M + x.y;
	ftemp = ( j * x.y - M );
	S = S + 1.0f / ( j * ( j - 1.0f ) ) * ftemp * ftemp;

	for (int i = 1; i < MSD_ELEM_PER_THREAD; i++)
	{
		pos = pos + WARP;
		x = __ldg(&d_input[pos]);

		j = j + 1.0f;
		M = M + x.x;
		ftemp = ( j * x.x - M );
		S = S + 1.0f / ( j * ( j - 1.0f ) ) * ftemp * ftemp;

		j = j + 1.0f;
		M = M + x.y;
		ftemp = ( j * x.y - M );
		S = S + 1.0f / ( j * ( j - 1.0f ) ) * ftemp * ftemp;
	}
	Ms[threadIdx.x] = M;
	Ss[threadIdx.x] = S;

	__syncthreads();

	// now all threads had saved their work, reduction follows

	// first we must load initial values
	//j=2*MSD_ELEM_PER_THREAD; // value of j is preserved during kernel's execution
	for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1)
	{
		if (threadIdx.x < i)
		{
			j = j * 2;
			ftemp = ( M - Ms[i + threadIdx.x] );
			S = S + Ss[i + threadIdx.x] + ( 1.0f / j ) * ftemp * ftemp;
			M = M + Ms[i + threadIdx.x];

			Ms[threadIdx.x] = M;
			Ss[threadIdx.x] = S;
		}
		__syncthreads();
	}

	// by now we should have only 32 partial results. shuffle reduction follows
	for (int q = HALF_WARP; q > 0; q = q >> 1)
	{
		j = j * 2;
		ftemp = ( M - __shfl_down(M, q) );
		S = S + __shfl_down(S, q) + ( 1.0f / j ) * ftemp * ftemp;
		M = M + __shfl_down(M, q);
	}

	//----------------------------------------------
	//---- Writing data
	if (warp_id == 0 && threadIdx.x == 0)
	{
		d_output[2 * blockIdx.x] = M;
		d_output[2 * blockIdx.x + 1] = S;
	}
}

__global__ void MSD_GPU_remainder(float const* __restrict__ d_input, float *d_output, int remainder)
{
	__shared__ float Ms[WARP * MSD_WARPS_PER_BLOCK];
	__shared__ float Ss[WARP * MSD_WARPS_PER_BLOCK];
	__shared__ float js[WARP * MSD_WARPS_PER_BLOCK];

	int warp_id, pos;
	float x;
	float M;
	float S;
	float j, jv;
	float ftemp;

	warp_id = threadIdx.x >> 5;

	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos = threadIdx.x;
	if (remainder > blockDim.x)
	{
		M = __ldg(&d_input[pos]);
		S = 0;
		j = 1.0f;
		pos = pos + blockDim.x;
		while (pos < remainder)
		{
			x = __ldg(&d_input[pos]);
			j = j + 1.0f;
			M = M + x;
			ftemp = ( j * x - M );
			S = S + 1.0f / ( j * ( j - 1.0f ) ) * ftemp * ftemp;
			pos = pos + blockDim.x;
		}

		Ms[threadIdx.x] = M;
		Ss[threadIdx.x] = S;
		js[threadIdx.x] = j;

		__syncthreads();

		// first we must load initial values
		for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1)
		{
			if (threadIdx.x < i)
			{
				jv = js[i + threadIdx.x];
				ftemp = ( jv / j * M - Ms[i + threadIdx.x] );
				S = S + Ss[i + threadIdx.x] + ( j / ( jv * ( j + jv ) ) ) * ftemp * ftemp;
				M = M + Ms[i + threadIdx.x];
				j = j + jv;

				Ms[threadIdx.x] = M;
				Ss[threadIdx.x] = S;
				js[threadIdx.x] = j;
			}
			__syncthreads();
		}

		// by now we should have only 32 partial results. shuffle reduction follows
		for (int q = HALF_WARP; q > 0; q = q >> 1)
		{
			jv = __shfl_down(j, q);
			ftemp = ( jv / j * M - __shfl_down(M, q) );
			S = S + __shfl_down(S, q) + ( j / ( jv * ( j + jv ) ) ) * ftemp * ftemp;
			M = M + __shfl_down(M, q);
			j = j + jv;
		}

	}
	else
	{
		if (threadIdx.x == 0)
		{	// This assumes remainder to be small < 32
			pos = 0;
			M = __ldg(&d_input[pos]);
			S = 0;
			j = 1.0f;
			for (pos = 1; pos < remainder; pos++)
			{
				x = __ldg(&d_input[pos]);
				j = j + 1.0f;
				M = M + x;
				ftemp = ( j * x - M );
				S = S + 1.0f / ( j * ( j - 1.0f ) ) * ftemp * ftemp;
			}
		}
	}

	//----------------------------------------------
	//---- Writing data
	if (warp_id == 0 && threadIdx.x == 0)
	{
		d_output[0] = M;
		d_output[1] = S;
		d_output[2] = j;
	}
}

__global__ void MSD_GPU_final(float *d_input, float *d_output, int size, int tail, float nElements)
{
	__shared__ float Ms[WARP * MSD_WARPS_PER_BLOCK];
	__shared__ float Ss[WARP * MSD_WARPS_PER_BLOCK];
	__shared__ float js[WARP * MSD_WARPS_PER_BLOCK];

	int warp_id, pos;
	float M;
	float S;
	float j, jv;
	float ftemp;

	warp_id = threadIdx.x >> 5;

	//----------------------------------------------
	//---- Calculating of streaming mean and sum of squares
	pos = threadIdx.x;
	if (size > blockDim.x)
	{
		M = d_input[2 * pos];
		S = d_input[2 * pos + 1];
		j = 2 * WARP * MSD_ELEM_PER_THREAD * MSD_WARPS_PER_BLOCK;
		jv = 2 * WARP * MSD_ELEM_PER_THREAD * MSD_WARPS_PER_BLOCK;
		pos = pos + blockDim.x;
		while (pos < size)
		{
			ftemp = ( jv / j * M - d_input[2 * pos] );
			S = S + d_input[2 * pos + 1] + ( j / ( jv * ( j + jv ) ) ) * ftemp * ftemp;
			M = M + d_input[2 * pos];
			j = j + jv;
			pos = pos + blockDim.x;
		}

		__syncthreads();

		Ms[threadIdx.x] = M;
		Ss[threadIdx.x] = S;
		js[threadIdx.x] = j;

		// now all threads had saved their work, reduction follows

		// first we must load initial values
		for (int i = ( blockDim.x >> 1 ); i > HALF_WARP; i = i >> 1)
		{
			if (threadIdx.x < i)
			{
				jv = js[i + threadIdx.x];
				ftemp = ( jv / j * M - Ms[i + threadIdx.x] );
				S = S + Ss[i + threadIdx.x] + ( j / ( jv * ( j + jv ) ) ) * ftemp * ftemp;
				M = M + Ms[i + threadIdx.x];
				j = j + jv;

				Ms[threadIdx.x] = M;
				Ss[threadIdx.x] = S;
				js[threadIdx.x] = j;
			}
			__syncthreads();
		}

		// by now we should have only 32 partial results. shuffle reduction follows
		for (int q = HALF_WARP; q > 0; q = q >> 1)
		{
			jv = __shfl_down(j, q);
			ftemp = ( jv / j * M - __shfl_down(M, q) );
			S = S + __shfl_down(S, q) + ( j / ( jv * ( j + jv ) ) ) * ftemp * ftemp;
			M = M + __shfl_down(M, q);
			j = j + jv;
		}

		if (tail > 0 && threadIdx.x == 0)
		{
			jv = d_input[2 * size + 2];
			ftemp = ( jv / j * M - d_input[2 * size] );
			S = S + d_input[2 * size + 1] + ( j / ( jv * ( j + jv ) ) ) * ftemp * ftemp;
			M = M + d_input[2 * size];
			j = j + jv;
		}
	}
	else
	{
		if (threadIdx.x == 0)
		{
			pos = 0;
			M = d_input[2 * pos];
			S = d_input[2 * pos + 1];
			j = 2 * WARP * MSD_ELEM_PER_THREAD * MSD_WARPS_PER_BLOCK;
			for (pos = 1; pos < size; pos++)
			{
				j = j * 2;
				ftemp = ( M - d_input[2 * pos] );
				S = S + d_input[2 * pos + 1] + ( 1.0f / j ) * ftemp * ftemp;
				M = M + d_input[2 * pos];
			}

			if (tail > 0)
			{
				jv = d_input[2 * size + 2];
				ftemp = ( jv / j * M - d_input[2 * size] );
				S = S + d_input[2 * size + 1] + ( j / ( jv * ( j + jv ) ) ) * ftemp * ftemp;
				M = M + d_input[2 * size];
				j = j + jv;
			}
		}
	}

	//----------------------------------------------
	//---- Writing data
	if (warp_id == 0 && threadIdx.x == 0)
	{
		d_output[0] = M / nElements;
		d_output[1] = sqrt(S / nElements);
		d_output[2] = nElements;

	}
}

#endif
