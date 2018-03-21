#ifndef ASTROACCELERATE_GPUANALYSIS_H_
#define ASTROACCELERATE_GPUANALYSIS_H_

//void analysis_GPU(float *h_peak_list, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int outBin, int *maxshift, int max_ndms, int *ndms, float cutoff, float sigma_constant, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, cudaStream_t streams, int candidate_algorithm, int enable_sps_baselinenoise, float *d_MSD_workarea, unsigned short *d_output_taps, float *d_MSD_interpolated, float *h_MSD_DIT, float *h_MSD_interpolated, int *gmem_peak_pos, unsigned long int maxTimeSamp);
void analysis_GPU(float *h_peak_list, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int outBin, int *maxshift, int max_ndms, int *ndms, float cutoff, float sigma_constant, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, cudaStream_t streams, int candidate_algorithm, int enable_sps_baselinenoise, float *d_MSD_workarea, unsigned short *d_output_taps, float *d_MSD_interpolated, int *gmem_peak_pos, unsigned long int maxTimeSamp);

#endif





