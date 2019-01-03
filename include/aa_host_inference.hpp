#ifndef ASTRO_ACCELERATE_AA_HOST_INFERENCE_HPP
#define ASTRO_ACCELERATE_AA_HOST_INFERENCE_HPP

namespace astroaccelerate {

  void cpu_blocked_bootstrap(float *data_array, unsigned int num_els, unsigned int num_bins, unsigned int num_boots, float *mean_data, float *cpu_mean, float *cpu_sd);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_HOST_INFERENCE_HPP

