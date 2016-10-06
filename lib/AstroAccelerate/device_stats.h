#ifndef __STATS_GPU__
#define __STATS_GPU__

extern void stats_gpu(cudaEvent_t event, cudaStream_t stream, int samps, float *mean, float *stddev, float *h_signal_power, float *d_signal_power);

#endif
