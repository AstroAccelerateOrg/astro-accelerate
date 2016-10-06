#ifndef __SINGLE_PULSE_SEARCH__
#define __SINGLE_PULSE_SEARCH__

extern void PD_SEARCH_init(void);
extern int PD_SEARCH(float *d_input, float *d_output, float *d_output_taps, float *d_MSD, int maxTaps, int nDMs, int nTimesamples);

#endif

