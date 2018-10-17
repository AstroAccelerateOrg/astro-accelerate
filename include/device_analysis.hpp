#ifndef ASTRO_ACCELERATE_DEVICE_ANALYSIS_HPP
#define ASTRO_ACCELERATE_DEVICE_ANALYSIS_HPP

#include "device_SPS_plan.hpp"

//void analysis_GPU(unsigned int *h_peak_list_DM, unsigned int *h_peak_list_TS, float *h_peak_list_SNR, unsigned int *h_peak_list_BW, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int outBin, int *maxshift, int max_ndms, int *ndms, float cutoff, float OR_sigma_multiplier, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, int candidate_algorithm, float *d_MSD_workarea, unsigned short *d_output_taps, float *d_MSD_interpolated, unsigned long int maxTimeSamples, int enable_sps_baselinenoise);
//void analysis_GPU(float *h_peak_list, size_t *peak_pos, size_t max_peak_size, SPS_DataDescription SPS_data, float *d_SPS_input, SPS_Parameters *SPS_params, MSD_Parameters *MSD_params);
void analysis_GPU(bool verbose, float *d_SPS_input, float *h_candidate_list, size_t &number_candidates, SPS_Plan spsplan);

#endif





