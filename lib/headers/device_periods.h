#ifndef ASTROACCELERATE_GPUPERIODICITY_SEARCH_H_
#define ASTROACCELERATE_GPUPERIODICITY_SEARCH_H_

extern void GPU_periodicity(int range, int nsamp, int max_ndms, int processed, float sigma_cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp, int nHarmonics, int candidate_algorithm, int enable_outlier_rejection, float bln_sigma_constant, float *h_MSD_DIT, float *h_MSD_interpolated, cudaStream_t streams);
//extern void GPU_periodicity(int range, int nsamp, int max_ndms, int processed, float sigma_cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp, int nHarmonics, int candidate_algorithm, int enable_outlier_rejection, float bln_sigma_constant, cudaStream_t streams);

#endif








