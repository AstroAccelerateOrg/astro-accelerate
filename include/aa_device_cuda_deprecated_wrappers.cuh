#ifndef ASTRO_ACCELERATE_AA_DEVICE_CUDA_DEPRECATED_WRAPPERS_CUH
#define ASTRO_ACCELERATE_AA_DEVICE_CUDA_DEPRECATED_WRAPPERS_CUH

namespace astroaccelerate {

#ifndef AA_ASSUME_MASK
#define AA_ASSUME_MASK 0xFFFFFFFF
#endif

  /** \brief Wrapper function that implements CUDA shfl for old and new CUDA implementation. */
  template<typename T, typename U, typename V>
  __device__ __inline__ T aa_shfl(const U &mask, const T &XX, const V &YY) {
#if(CUDART_VERSION >= 9000)
    return(__shfl_sync(mask, XX, YY));
#else
    return(__shfl(XX, YY));
#endif
  }

  /** \brief Wrapper function that implements CUDA shfl_up for old and new CUDA implementation. */
  template<typename T, typename U, typename V>
  __device__ __inline__ T aa_shfl_up(const U &mask, const T &XX, const V &YY) {
#if(CUDART_VERSION >= 9000)
    return(__shfl_up_sync(mask, XX, YY));
#else
    return(__shfl_up(XX, YY));
#endif
  }

  /** \brief Wrapper function that implements CUDA shfl_down for old and new CUDA implementation. */
  template<typename T, typename U, typename V>
  __device__ __inline__ T aa_shfl_down(const U &mask, const T &XX, const V &YY) {
#if(CUDART_VERSION >= 9000)
    return(__shfl_down_sync(mask, XX, YY));
#else
    return(__shfl_down(XX, YY));
#endif
  }

  /** \brief Wrapper function that implements CUDA shfl_xor for old and new CUDA implementation. */
  template<typename T, typename U, typename V>
  __device__ __inline__ T aa_shfl_xor(const U &mask, const T &XX, const V &YY) {
#if(CUDART_VERSION >= 9000)
    return(__shfl_xor_sync(mask, XX, YY));
#else
    return(__shfl_xor(XX, YY));
#endif
  }

  /** \brief Wrapper function that implements CUDA ballot for old and new CUDA implementation. */
  template<typename T, typename U>
  __device__ __inline__ T aa_ballot(const U &mask, const T &XX) {
#if(CUDART_VERSION >= 9000)
    return(__ballot_sync(mask, XX));
#else
    return(__ballot(XX));
#endif
  }

} //namespace astroaccelerate

#endif //ASTRO_ACCELERATE_AA_DEVICE_CUDA_DEPRECATED_WRAPPERS_CUH
