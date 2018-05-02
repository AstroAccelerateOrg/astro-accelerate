#ifndef ASTROACCELERATE_GPUANALYSIS_H_
#define ASTROACCELERATE_GPUANALYSIS_H_

//void analysis_GPU(float *h_peak_list, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int outBin, int *maxshift, int max_ndms, int *ndms, float cutoff, float sigma_constant, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, int candidate_algorithm, int enable_sps_baselinenoise);

void analysis_GPU( float *h_peak_list, size_t *peak_pos, size_t max_peak_size, SPS_Data_Description SPS_data, float *d_SPS_input, SPS_Parameters *SPS_params, MSD_Parameters *MSD_params, AA_Parameters *AA_params);

#endif





