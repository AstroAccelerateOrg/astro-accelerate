#ifndef ASTROACCELERATE_POWER_GPU_H_
#define ASTROACCELERATE_POWER_GPU_H_

#include <cufft.h>

extern void power_gpu(cudaEvent_t event, cudaStream_t stream, int samps, int acc, cufftComplex *d_signal_fft, float *d_signal_power);

#endif

