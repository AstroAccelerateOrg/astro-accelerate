// James Sharpe's peak finding code

#ifndef __PEAK_FIND__
#define __PEAK_FIND__

#include<vector>
#include "device_SPS_BC_plan.h"

extern void PEAK_FIND(float *d_output_SNR, ushort *d_output_taps, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration);
extern void PEAK_FIND_FOR_FDAS(float *d_ffdot_plane, float *d_peak_list, float *d_MSD, int nDMs, int nTimesamples, float threshold, unsigned int max_peak_size, unsigned int *gmem_peak_pos, float DM_trial);

extern void Peak_find_for_periodicity_search_old(float *d_input_SNR, ushort *d_input_harmonics, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, float *d_MSD, int DM_shift, int inBin);

extern void Peak_find_for_periodicity_search(float *d_input_SNR, ushort *d_input_harmonics, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, float const* __restrict__ d_MSD, int DM_shift, int inBin);


#endif
