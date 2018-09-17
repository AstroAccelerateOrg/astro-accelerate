#ifndef ASTRO_ACCELERATE_DEVICE_DEDISPERSE_HPP
#define ASTRO_ACCELERATE_DEVICE_DEDISPERSE_HPP

extern void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, unsigned short *d_input, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms, int nbits, int failsafe);

#endif

