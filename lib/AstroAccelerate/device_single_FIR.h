#ifndef ASTROACCELERATE_SINGLE_FIR_H_
#define ASTROACCELERATE_SINGLE_FIR_H_

extern void PD_FIR_init(void);
extern int PD_FIR(float *d_input, float *d_output, int nTaps, int nDMs, int nTimesamples);

#endif
