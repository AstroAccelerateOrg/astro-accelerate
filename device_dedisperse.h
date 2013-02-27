#ifndef __DEDISPERSION__
#define __DEDISPERSION__

extern void dedisperse(size_t inputsize, float *d_input, size_t outputsize, float *d_output, int nchans, int nsamp, int maxshift, float dm_low, int ndms, int kernel_type, float tsamp, float dm_step);

#endif

