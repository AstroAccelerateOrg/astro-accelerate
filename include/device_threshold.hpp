#ifndef ASTRO_ACCELERATE_DEVICE_THRESHOLD_HPP
#define ASTRO_ACCELERATE_DEVICE_THRESHOLD_HPP

#include <tuple>
#include "device_BC_plan.hpp"
#include "device_SPS_plan.hpp"

extern void THR_init(void);
extern int SPDT_threshold(float *d_input, unsigned short *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nDMs, int nTimesamples, int shift, SPS_Plan spsplan, int max_iteration, int max_list_size, std::tuple<float, float, float> &dmlimits, float sampling_time, float inBin, float start_time);

extern int Threshold_for_periodicity_old(float *d_input, unsigned short *d_input_harms, float *d_output_list, int *gmem_pos, float *d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int inBin, int max_list_size);

extern int Threshold_for_periodicity(float *d_input, unsigned short *d_input_harms, float *d_output_list, int *gmem_pos, float *d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int inBin, int max_list_size);

#endif

