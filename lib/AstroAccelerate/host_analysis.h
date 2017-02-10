#ifndef SKA_ASTROACCELERATE_ANALYSIS_H_
#define SKA_ASTROACCELERATE_ANALYSIS_H_

#include <vector>

void analysis(int i, float tstart, int t_processed, int nsamp, int nchans, int maxshift, int max_ndms, int *ndms, int *outBin, float cutoff, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, std::vector<float> &output_sps, int analysis_call);

#endif

