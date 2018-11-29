/**
 * Example code for linking against astro-accelerate library.
 *
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion.cpp -o test
 */

#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_permitted_pipelines_2.hpp"

#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"

using namespace astroaccelerate;

int main() {
  aa_ddtr_plan ddtr_plan;
  ddtr_plan.add_dm(0, 370, 0.307, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).
  ddtr_plan.add_dm(370, 740, 0.652, 2, 2);

  // Filterbank metadata
  // (Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs")
  const double tstart = 50000;
  const double tsamp = 6.4E-5;
  const double nbits = 8;
  const double nsamples = 937984;
  const double fch1 = 1564;
  const double foff = -0.208984;
  const double nchans = 2048;
  
  aa_filterbank_metadata metadata(tstart, tsamp, nbits, nsamples, fch1, foff, nchans);
  
  const size_t free_memory = 2147483648; // Free memory on the GPU in bytes
  bool enable_analysis = true;       // The strategy will be optimised to run just dedispersion
  aa_ddtr_strategy ddtr_strategy(ddtr_plan, metadata, free_memory, enable_analysis);
  
  if(!(ddtr_strategy.ready())) {
    std::cout << "ERROR: ddtr_strategy not ready." << std::endl;
    return 0;
  }

  std::vector<unsigned short> input_data(nsamples*nchans);

  for(auto& i : input_data) {
    i = 0.0;
  }

  const float sigma_cutoff = 6.0;
  const float sigma_constant = 4.0;
  const float max_boxcar_width_in_sec = 0.05;
  const aa_analysis_plan::selectable_candidate_algorithm algo = aa_analysis_plan::selectable_candidate_algorithm::off;
  
  aa_analysis_plan analysis_plan(ddtr_strategy, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, algo, false);
  aa_analysis_strategy analysis_strategy(analysis_plan);

  if(!(analysis_strategy.ready())) {
    std::cout << "ERROR: analysis_strategy not ready." << std::endl;
    return 0;
  }
  
  aa_permitted_pipelines_2<aa_compute::module_option::zero_dm, false> runner(ddtr_strategy, analysis_strategy, input_data.data());
  if(runner.setup()) {
    while(runner.next()) {
      std::cout << "NOTICE: Pipeline running over next chunk." << std::endl;
    }
  }
  
  std::cout << "NOTICE: Finished." << std::endl;
  return 0;
}
