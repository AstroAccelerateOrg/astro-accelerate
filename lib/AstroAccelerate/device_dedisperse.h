#ifndef __DEDISPERSION__
#define __DEDISPERSION__

extern void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, float *d_input, cudaTextureObject_t tex, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms);

#endif

