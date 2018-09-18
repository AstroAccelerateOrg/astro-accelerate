#ifndef ASTRO_ACCELERATE_DEVICE_STRETCH_HPP
#define ASTRO_ACCELERATE_DEVICE_STRETCH_HPP

extern void stretch_gpu(cudaEvent_t event, cudaStream_t stream, int acc, int samps, float tsamp, float *d_input, float *d_output);

#endif

