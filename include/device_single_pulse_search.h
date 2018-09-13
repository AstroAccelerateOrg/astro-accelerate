#ifndef ASTROACCELERATE_SINGLE_PULSE_SEARCH_H_
#define ASTROACCELERATE_SINGLE_PULSE_SEARCH_H_

extern void PD_SEARCH_init(void);
extern int PD_SEARCH(float *d_input, float *d_output, float *d_output_taps, float *d_MSD, int maxTaps, int nDMs, int nTimesamples);

#endif

