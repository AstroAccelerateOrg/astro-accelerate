#ifndef ASTROACCELERATE_THRESHOLD_H_
#define ASTROACCELERATE_THRESHOLD_H_

#include<vector>
#include "device_SPS_PD_plan.h"


extern void THR_init(void);
extern int THRESHOLD(float *d_input, ushort *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nDMs, int nTimesamples, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int max_list_size, float dm_step, float dm_low, float sampling_time, float inBin, float start_time);

extern int Threshold_for_periodicity_old(float *d_input, ushort *d_input_harms, float *d_output_list, int *gmem_pos, float *d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int inBin, int max_list_size);

extern int Threshold_for_periodicity(float *d_input, ushort *d_input_harms, float *d_output_list, int *gmem_pos, float const* __restrict__ d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int inBin, int max_list_size);

#endif

