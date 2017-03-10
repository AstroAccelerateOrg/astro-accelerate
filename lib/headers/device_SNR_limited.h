#ifndef ASTROACCELERATE_SNR_LIMITED_H_
#define ASTROACCELERATE_SNR_LIMITED_H_

extern void SNR_limited_init(void);
extern int SNR_limited(float *d_FIR_input, float *d_SNR_output, float *d_SNR_taps, float *d_MSD, int nTaps, int nDMs, int nTimesamples, int offset);

#endif

