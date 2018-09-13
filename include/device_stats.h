#ifndef ASTROACCELERATE_STATS_GPU_H_
#define ASTROACCELERATE_STATS_GPU_H_

extern void stats_gpu(cudaEvent_t event, cudaStream_t stream, int samps, float *mean, float *stddev, float *h_signal_power, float *d_signal_power);

#endif
