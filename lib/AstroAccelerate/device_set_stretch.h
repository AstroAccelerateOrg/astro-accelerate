#ifndef __SET_STRETCH_GPU__
#define __SET_STRETCH_GPU__

extern void set_stretch_gpu(cudaEvent_t event, cudaStream_t stream, int samps, float mean, float *d_input);

#endif

