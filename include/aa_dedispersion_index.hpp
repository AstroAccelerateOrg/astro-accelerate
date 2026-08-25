#ifndef ASTRO_ACCELERATE_AA_DEDISPERSION_INDEX_HPP
#define ASTRO_ACCELERATE_AA_DEDISPERSION_INDEX_HPP

#include <cstddef>

#if defined(__CUDACC__)
#define AA_DEDISPERSION_HOST_DEVICE __host__ __device__
#else
#define AA_DEDISPERSION_HOST_DEVICE
#endif

namespace astroaccelerate {

AA_DEDISPERSION_HOST_DEVICE constexpr std::size_t
dedispersion_flattened_index(std::size_t outer_index,
                             std::size_t stride,
                             std::size_t inner_index) {
  return outer_index * stride + inner_index;
}

} // namespace astroaccelerate

#undef AA_DEDISPERSION_HOST_DEVICE

#endif // ASTRO_ACCELERATE_AA_DEDISPERSION_INDEX_HPP
