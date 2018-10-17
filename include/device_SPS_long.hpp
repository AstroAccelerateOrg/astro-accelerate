#ifndef ASTRO_ACCELERATE_DEVICE_SPS_LONG_HPP
#define ASTRO_ACCELERATE_DEVICE_SPS_LONG_HPP

#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include <helper_cuda.h>

#include "device_BC_plan.hpp"
#include "params.hpp"
#include "device_BC_plan.hpp"
#include "device_SPS_long_kernel.hpp"
#include "device_SPS_plan.hpp"

size_t Get_memory_requirement_of_SPS();
extern void PD_SEARCH_LONG_init();
extern int SPDT_search_long_MSD_plane(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD_interpolated, SPS_Plan spsplan, int max_iteration, int nTimesamples, int nDMs);	

#endif
