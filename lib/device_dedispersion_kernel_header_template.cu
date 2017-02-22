#ifndef DEDISPERSE_KERNEL_H_
#define DEDISPERSE_KERNEL_H_

#define ARRAYSIZE SDIVINT * SDIVINDM

#include "float.h"
#include "AstroAccelerate/kernel_params.h"

// Stores temporary shift values
__device__ __constant__ float dm_shifts[8192];
__device__ __constant__ int i_nsamp, i_nchans, i_t_processed_s;
//__device__  __shared__ ushort2 f_line[UNROLLS][ARRAYSIZE + 1];

//{{{ shared_dedisperse_loop

