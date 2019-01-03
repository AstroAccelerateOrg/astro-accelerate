#ifndef ASTRO_ACCELERATE_AA_DEVICE_STATS_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_STATS_HPP

namespace astroaccelerate {

  extern void stats_gpu(cudaEvent_t event, cudaStream_t stream, int samps, float *mean, float *stddev, float *h_signal_power, float *d_signal_power);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_STATS_HPP
