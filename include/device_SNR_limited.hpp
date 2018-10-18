#ifndef ASTRO_ACCELERATE_DEVICE_SNR_LIMITED_HPP
#define ASTRO_ACCELERATE_DEVICE_SNR_LIMITED_HPP

#include "params.hpp"
#include "device_SNR_limited_kernel.hpp"

extern void SNR_limited_init(void);
extern int SNR_limited(float *d_FIR_input, float *d_SNR_output, float *d_SNR_taps, float *d_MSD, int nTaps, int nDMs, int nTimesamples, int offset);

#endif

