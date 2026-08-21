#include <cstddef>
#include <iostream>


#include "aa_pulscan_runner.hpp"

using namespace astroaccelerate;

int main() {
  std::cout << "Running test_pulscan_compatibility_helpers.cpp" << std::endl;

  struct test_case {
    size_t input;
    size_t expected;
  };

  const test_case cases[] = {
      {0, 0},
      {4095, 0},
      {4096, 4096},
      {5000, 4096},
      {8192, 8192},
      {8193, 8192},
      {12288, 12288},
  };

  for(const auto& test : cases) {
    const size_t actual = pulscan_truncate_fft_bins(test.input);
    if(actual != test.expected) {
      std::cout << "truncate(" << test.input << ") expected " << test.expected
                << " but got " << actual << std::endl;
      return 1;
    }
  }

  struct series_length_case {
    size_t available_samples;
    size_t planned_samples;
    size_t expected;
  };

  const series_length_case series_length_cases[] = {
      {0, 0, 0},
      {0, 4096, 0},
      {4096, 0, 4096},
      {4096, 8192, 4096},
      {4096, 4096, 4096},
      {4096, 2048, 2048},
      {28256704, 16777216, 16777216},
  };

  for(const auto& test : series_length_cases) {
    const size_t actual = pulscan_select_samples_per_series(
        test.available_samples, test.planned_samples);
    if(actual != test.expected) {
      std::cout << "select_samples_per_series(" << test.available_samples
                << ", " << test.planned_samples << ") expected "
                << test.expected << " but got " << actual << std::endl;
      return 1;
    }
  }

  std::cout << "Runs" << std::endl;
  return 0;
}
