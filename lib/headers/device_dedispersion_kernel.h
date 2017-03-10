#ifndef ASTROACCELERATE_DEDISPERSION_KERNEL_H_
#define ASTROACCELERATE_DEDISPERSION_KERNEL_H_


//#define ARRAYSIZE DIVINT * DIVINDM
#define ARRAYSIZE SDIVINT * SDIVINDM

//#include <cuda.h>
//#include <cuda_runtime.h>
#include "params.h"

// Stores temporary shift values
__device__ __constant__ float dm_shifts[15500];
__device__ __constant__ int   i_nsamp, i_nchans, i_t_processed_s, i_t_processed_c;
//__device__ __shared__ float fa_line[ARRAYSIZE+1];
__device__ __shared__ float2 fa_line[ARRAYSIZE+1];
__device__ __shared__ float2 fb_line[ARRAYSIZE+1];
__device__ __shared__ float2 fc_line[ARRAYSIZE+1];
__device__ __shared__ float2 fd_line[ARRAYSIZE+1];
//texture<float,1,cudaReadModeElementType>inTex;


#endif

