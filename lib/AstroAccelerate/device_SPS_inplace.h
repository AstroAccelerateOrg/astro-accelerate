#ifndef __SPS_INPLACE__
#define __SPS_INPLACE__

extern void PD_SEARCH_INPLACE_init(void);
extern int PD_SEARCH_INPLACE(float *d_input, unsigned char *d_output_taps, float *d_MSD, int maxTaps, int nDMs, int nTimesamples);

#endif

