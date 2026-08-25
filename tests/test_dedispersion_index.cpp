#include <climits>
#include <cstddef>
#include <iostream>

#include "aa_dedispersion_index.hpp"

using namespace astroaccelerate;

int main() {
  std::cout << "Running test_dedispersion_index.cpp" << std::endl;

  constexpr std::size_t samples_per_channel = 28256704;
  constexpr std::size_t last_safe_int_index =
    dedispersion_flattened_index(75, samples_per_channel, 0);
  constexpr std::size_t first_overflowing_int_index =
    dedispersion_flattened_index(76, samples_per_channel, 0);
  constexpr std::size_t last_input_index =
    dedispersion_flattened_index(319, samples_per_channel,
                                 samples_per_channel - 1);

  static_assert(last_safe_int_index <= INT_MAX,
                "The regression boundary should fit in a signed int");
  static_assert(first_overflowing_int_index > INT_MAX,
                "The regression boundary should require a 64-bit index");
  static_assert(first_overflowing_int_index == 2147509504ULL,
                "The flattened channel offset must not wrap");
  static_assert(last_input_index > UINT_MAX,
                "Production-sized input indexing must exceed 32 bits");
  static_assert(last_input_index == 9042145279ULL,
                "The full production-sized input offset must not truncate");

  const std::size_t advanced = dedispersion_flattened_index(
    4, samples_per_channel, last_safe_int_index);
  if(advanced !=
     dedispersion_flattened_index(79, samples_per_channel, 0)) {
    std::cout << "Channel-stride accumulation wrapped or truncated"
              << std::endl;
    return 1;
  }

  std::cout << "Runs" << std::endl;
  return 0;
}
