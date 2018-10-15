#ifndef ASTRO_ACCELERATE_DEVICE_STATS_HPP
#define ASTRO_ACCELERATE_DEVICE_STATS_HPP

extern void stats_gpu(cudaEvent_t  event,
                      cudaStream_t stream,
                      int          samps,
                      float *      mean,
                      float *      stddev,
                      float *      h_signal_power,
                      float *      d_signal_power);

#endif
