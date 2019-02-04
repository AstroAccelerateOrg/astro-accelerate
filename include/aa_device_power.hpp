#ifndef ASTRO_ACCELERATE_AA_DEVICE_POWER_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_POWER_HPP

#include <cufft.h>

namespace astroaccelerate {

extern void power_gpu(cudaEvent_t event, cudaStream_t stream, int samps, int acc, cufftComplex *d_signal_fft, float *d_signal_power);
extern void power_and_interbin_gpu(float2 *d_input_complex, float *d_output_power, float *d_output_interbinning, int nTimesamples, int nDMs);
extern void simple_power_and_interbin(float2 *d_input_complex, float *d_output_power, float *d_output_interbinning, int nTimesamples, int nDMs);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_POWER_HPP

