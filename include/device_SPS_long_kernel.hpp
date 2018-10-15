#ifndef ASTRO_ACCELERATE_DEVICE_SPS_LONG_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_SPS_LONG_KERNEL_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include "device_SPS_inplace_kernel.hpp"
#include "params.hpp"

void call_kernel_SPDT_GPU_1st_plane(dim3          grid_size,
                                    dim3          block_size,
                                    float const * d_input,
                                    float *       d_bv_out,
                                    float *       d_decimated,
                                    float *       d_output_SNR,
                                    ushort *      d_output_taps,
                                    float2 const *d_MSD,
                                    int           nTimesamples,
                                    int           nBoxcars,
                                    const int     dtm);

void call_kernel_SPDT_GPU_Nth_plane(dim3          grid_size,
                                    dim3          block_size,
                                    float const * d_input,
                                    float *       d_bv_in,
                                    float *       d_bv_out,
                                    float *       d_decimated,
                                    float *       d_output_SNR,
                                    ushort *      d_output_taps,
                                    float2 const *d_MSD,
                                    const int     nTimesamples,
                                    const int     nBoxcars,
                                    const int     startTaps,
                                    const int     DIT_value,
                                    const int     dtm);

#endif
