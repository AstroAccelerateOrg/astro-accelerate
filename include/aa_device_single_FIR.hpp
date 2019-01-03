#ifndef ASTRO_ACCELERATE_AA_DEVICE_SINGLE_FIR_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SINGLE_FIR_HPP

namespace astroaccelerate {

  extern void PD_FIR_init(void);
  extern int PD_FIR(float *d_input, float *d_output, int nTaps, int nDMs, int nTimesamples);
  extern int GPU_FIRv1_wrapper(float *d_input, float *d_output, int nTaps, unsigned int nDMs, unsigned int nTimesamples);
  extern int PPF_L1(float *d_input, float *d_output, int nChannels, int nSpectra, int nTaps);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SINGLE_FIR_HPP
