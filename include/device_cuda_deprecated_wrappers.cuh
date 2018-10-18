#ifndef ASTRO_ACCELERATE_DEVICE_CUDA_DEPRECATED_WRAPPERS_CUH
#define ASTRO_ACCELERATE_DEVICE_CUDA_DEPRECATED_WRAPPERS_CUH

#ifndef AA_ASSUME_MASK
#define AA_ASSUME_MASK 0xFFFFFFFF
#endif

template<typename T, typename U, typename V>
__device__ __inline__ T aa_shfl(const U &mask, const T &XX, const V &YY) {
#if(CUDART_VERSION >= 9000)
  return(__shfl_sync(mask, XX, YY));
#else
  return(__shfl(XX, YY));
#endif
}

template<typename T, typename U, typename V>
__device__ __inline__ T aa_shfl_up(const U &mask, const T &XX, const V &YY) {
#if(CUDART_VERSION >= 9000)
  return(__shfl_up_sync(mask, XX, YY));
#else
  return(__shfl_up(XX, YY));
#endif
}

template<typename T, typename U, typename V>
__device__ __inline__ T aa_shfl_down(const U &mask, const T &XX, const V &YY) {
#if(CUDART_VERSION >= 9000)
  return(__shfl_down_sync(mask, XX, YY));
#else
  return(__shfl_down(XX, YY));
#endif
}

template<typename T, typename U, typename V>
__device__ __inline__ T aa_shfl_xor(const U &mask, const T &XX, const V &YY) {
#if(CUDART_VERSION >= 9000)
  return(__shfl_xor_sync(mask, XX, YY));
#else
  return(__shfl_xor(XX, YY));
#endif
}

template<typename T, typename U>
__device__ __inline__ T aa_ballot(const U &mask, const T &XX) {
#if(CUDART_VERSION >= 9000)
  return(__ballot_sync(mask, XX));
#else
  return(__ballot(XX));
#endif
}

#endif //ASTRO_ACCELERATE_DEVICE_CUDA_DEPRECATED_WRAPPERS_CUH
