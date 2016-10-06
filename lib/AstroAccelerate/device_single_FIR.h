#ifndef __SINGLE_FIR__
#define __SINGLE_FIR__

extern void PD_FIR_init(void);
extern int PD_FIR(float *d_input, float *d_output, int nTaps, int nDMs, int nTimesamples);

#endif
