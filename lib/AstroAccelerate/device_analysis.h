#ifndef ASTROACCELERATE_GPUANALYSIS_H_
#define ASTROACCELERATE_GPUANALYSIS_H_

void analysis_GPU(int i, float tstart, int t_processed, int nsamp, int nchans, int maxshift, int max_ndms, int *ndms, int *outBin, float cutoff, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp);

#endif





