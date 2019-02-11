#include <iostream>
#include "aa_periodicity_strategy.hpp"

using namespace astroaccelerate;

int main() {
  std::cout << "Running test_ddtr_strategy_1.cpp" << std::endl;
  const float  sigma_cutoff             = 0.0;
  const float  sigma_constant           = 0.0;
  const int    nHarmonics               = 0.0;
  const int    export_powers            = 0.0;
  const bool   candidate_algorithm      = false;
  const bool   enable_outlier_rejection = false;
  aa_periodicity_plan plan(sigma_cutoff,
			   sigma_constant,
			   nHarmonics,
			   export_powers,
			   candidate_algorithm,
			   enable_outlier_rejection);
  aa_periodicity_strategy strategy(plan);

  if(strategy.ready()) {
    // The above parameters should not allow strategy to become ready.
    std::cout << "Fail." << std::endl;
    return 0;
  }
  
  
  std::cout << "Runs" << std::endl;
  return 0;
}
