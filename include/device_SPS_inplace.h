#ifndef ASTROACCELERATE_SPS_INPLACE_H_
#define ASTROACCELERATE_SPS_INPLACE_H_

extern void PD_SEARCH_INPLACE_init(void);
extern int PD_SEARCH_INPLACE(float *d_input, unsigned char *d_output_taps, float *d_MSD, int maxTaps, int nDMs, int nTimesamples);

#endif

