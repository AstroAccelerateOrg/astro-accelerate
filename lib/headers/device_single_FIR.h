#ifndef ASTROACCELERATE_SINGLE_FIR_H_
#define ASTROACCELERATE_SINGLE_FIR_H_

extern void PD_FIR_init(void);
extern int PD_FIR(float *d_input, float *d_output, int nTaps, int nDMs, int nTimesamples);
extern int GPU_FIRv1_wrapper(float *d_input, float *d_output, int nTaps, unsigned int nDMs, unsigned int nTimesamples);
extern int PPF_L1(float *d_input, float *d_output, int nChannels, int nSpectra, int nTaps);

#endif
