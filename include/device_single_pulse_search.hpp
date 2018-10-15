#ifndef ASTRO_ACCELERATE_DEVICE_SINGLE_PULSE_SEARCH_HPP
#define ASTRO_ACCELERATE_DEVICE_SINGLE_PULSE_SEARCH_HPP

extern void PD_SEARCH_init(void);
extern int  PD_SEARCH(float *d_input,
                      float *d_output,
                      float *d_output_taps,
                      float *d_MSD,
                      int    maxTaps,
                      int    nDMs,
                      int    nTimesamples);

#endif
