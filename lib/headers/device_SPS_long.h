#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "device_SPS_PD_plan.h"

#ifndef __SPS_LONG__
#define __SPS_LONG__

size_t Get_memory_requirement_of_SPS();
extern void PD_SEARCH_LONG_init();
extern int SPDT_search_long_MSD_plane(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD_interpolated, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nTimesamples, int nDMs);	

#endif
