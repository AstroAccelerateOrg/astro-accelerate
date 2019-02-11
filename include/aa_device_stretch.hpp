#ifndef ASTRO_ACCELERATE_AA_DEVICE_STRETCH_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_STRETCH_HPP

namespace astroaccelerate {

  extern void stretch_gpu(cudaEvent_t event, cudaStream_t stream, int acc, int samps, float tsamp, float *d_input, float *d_output);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_STRETCH_HPP

