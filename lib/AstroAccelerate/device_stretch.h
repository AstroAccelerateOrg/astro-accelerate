#ifndef __STRETCH_GPU__
#define __STRETCH_GPU__

extern void stretch_gpu(cudaEvent_t event, cudaStream_t stream, int acc, int samps, float tsamp, float *d_input, float *d_output);

#endif

