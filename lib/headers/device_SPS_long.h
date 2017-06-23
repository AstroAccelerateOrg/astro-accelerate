#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "device_BC_plan.h"

#ifndef __SPS_LONG__
#define __SPS_LONG__

size_t Get_memory_requirement_of_SPS();
extern void PD_SEARCH_LONG_init();
extern int PD_SEARCH_LONG_BLN(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples);
extern int PD_SEARCH_LONG_BLN_EACH(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples, float sigma_constant);
extern int PD_SEARCH_LONG_LINAPPROX(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples);
extern int PD_SEARCH_LONG_LINAPPROX_EACH(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples);
int PD_SEARCH_LONG_BLN_LINAPPROX_EACH(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples, float sigma_constant);
#endif
