#ifndef __DEDISPERSION__
#define __DEDISPERSION__

<<<<<<< HEAD
extern void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, unsigned short *d_input, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms);
=======
extern void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, float *d_input, cudaTextureObject_t tex, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms);
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229

#endif

