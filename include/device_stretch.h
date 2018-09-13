#ifndef ASTROACCELERATE_STRETCH_GPU_H_
#define ASTROACCELERATE_STRETCH_GPU_H_

extern void stretch_gpu(cudaEvent_t event, cudaStream_t stream, int acc, int samps, float tsamp, float *d_input, float *d_output);

#endif

