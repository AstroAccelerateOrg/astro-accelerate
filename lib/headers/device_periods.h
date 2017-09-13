#ifndef ASTROACCELERATE_GPUPERIODICITY_SEARCH_H_
#define ASTROACCELERATE_GPUPERIODICITY_SEARCH_H_

extern void GPU_periodicity(int range, int nsamp, int max_ndms, int processed, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp, int nHarmonics);

#endif








