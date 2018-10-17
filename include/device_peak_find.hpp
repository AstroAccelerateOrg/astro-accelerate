// James Sharpe's peak finding code

#ifndef ASTRO_ACCELERATE_DEVICE_PEAK_FIND_HPP
#define ASTRO_ACCELERATE_DEVICE_PEAK_FIND_HPP

#include <tuple>
#include <npp.h>
#include <helper_cuda.h>

#include "params.hpp"
#include "device_peak_find_kernel.hpp"
#include "device_BC_plan.hpp"
#include "device_SPS_plan.hpp"

extern void SPDT_peak_find(float *d_output_SNR, unsigned short *d_output_taps, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, int shift, SPS_Plan spsplan, int max_iteration, std::tuple<float, float, float> dmlimits, float sampling_time, float inBin, float start_time);

extern void PEAK_FIND_FOR_FDAS(float *d_ffdot_plane, float *d_peak_list, float *d_MSD, int nDMs, int nTimesamples, float threshold, unsigned int max_peak_size, unsigned int *gmem_peak_pos, float DM_trial);

extern void Peak_find_for_periodicity_search_old(float *d_input_SNR, ushort *d_input_harmonics, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, float *d_MSD, int DM_shift, int inBin);

extern void Peak_find_for_periodicity_search(float *d_input_SNR, ushort *d_input_harmonics, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, float* d_MSD, int DM_shift, int inBin);


#endif
