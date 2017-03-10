#ifndef ASTROACCELERATE_DEBUG_H_
#define ASTROACCELERATE_DEBUG_H_

#include <time.h>

void debug(int test, clock_t start_time, int range, int *outBin, int enable_debug, int analysis, int output_dmt, int multi_file, float sigma_cutoff, float power, int max_ndms, float *user_dm_low, float *user_dm_high, float *user_dm_step, float *dm_low, float *dm_high, float *dm_step, int *ndms, int nchans, int nsamples, int nifs, int nbits, float tsamp, float tstart, float fch1, float foff, int maxshift, float max_dm, int nsamp, size_t gpu_inputsize, size_t gpu_outputsize, size_t inputsize, size_t outputsize);

#endif
