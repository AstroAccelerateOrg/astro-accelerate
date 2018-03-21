#ifndef ASTROACCELERATE_BIN_H_
#define ASTROACCELERATE_BIN_H_

extern void bin_gpu(unsigned short *d_input, float *d_output, int nchans, int nsamp, cudaStream_t streams);
extern int GPU_DiT_v2_wrapper(float *d_input, float *d_output, int nDMs, int nTimesamples, cudaStream_t streams);

#endif

