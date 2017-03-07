#ifndef ASTROACCELERATE_SET_STRETCH_GPU_H_
#define ASTROACCELERATE_SET_STRETCH_GPU_H_

extern void set_stretch_gpu(cudaEvent_t event, cudaStream_t stream, int samps, float mean, float *d_input);

#endif

