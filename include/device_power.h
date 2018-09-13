#ifndef ASTROACCELERATE_POWER_GPU_H_
#define ASTROACCELERATE_POWER_GPU_H_

#include <cufft.h>

extern void power_gpu(cudaEvent_t event, cudaStream_t stream, int samps, int acc, cufftComplex *d_signal_fft, float *d_signal_power);
extern void power_and_interbin_gpu(float2 *d_input_complex, float *d_output_power, float *d_output_interbinning, int nTimesamples, int nDMs);
extern void simple_power_and_interbin(float *d_input_complex, float *d_output_power, float *d_output_interbinning, int nTimesamples, int nDMs);

#endif

