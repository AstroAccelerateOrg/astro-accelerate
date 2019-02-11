#ifndef ASTRO_ACCELERATE_AA_DEVICE_THRESHOLD_SHARED_KERNEL_FUNCTIONS_CUH
#define ASTRO_ACCELERATE_AA_DEVICE_THRESHOLD_SHARED_KERNEL_FUNCTIONS_CUH

namespace astroaccelerate {

__device__ __inline__ float white_noise(float *value, float *hrms, float *mean, float *sd){
  return( __frsqrt_rn((*hrms)+1.0f)*((*value) - (*mean)*((*hrms)+1.0f))/((*sd)) );
}

__device__ __inline__ float inverse_white_noise(float *SNR, float *hrms, float *mean, float *sd){
  return( (*SNR)*__fsqrt_rn((*hrms)+1.0f)*(*sd) + (*mean)*((*hrms)+1.0f) );
}

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_THRESHOLD_SHARED_KERNEL_FUNCTIONS_CUH
