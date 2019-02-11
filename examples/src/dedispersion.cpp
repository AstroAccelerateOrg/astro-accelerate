/**
 * Example code for linking against astro-accelerate library.
 *
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion.cpp -o test
 */

#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_permitted_pipelines_1.hpp"

using namespace astroaccelerate;

int main() {
  aa_ddtr_plan ddtr_plan;
  ddtr_plan.add_dm(0, 370, 0.307, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).
  ddtr_plan.add_dm(370, 740, 0.652, 1, 1);
  ddtr_plan.add_dm(740, 1480, 1.266, 4, 4);
  ddtr_plan.add_dm(1480, 2950, 25.12, 8, 8);
  ddtr_plan.add_dm(2950, 5000, 4.000, 16, 16);

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
  bool enable_analysis = false;       // The strategy will be optimised to run just dedispersion
  aa_ddtr_strategy strategy(ddtr_plan, metadata, free_memory, enable_analysis);
  
  if(!(strategy.ready())) {
    std::cout << "There was an error" << std::endl;
    return 0;
  }

  std::vector<unsigned short> input_data(nsamples*nchans);

  for(auto &i : input_data) {
    i = 0.0;
  }
  
  aa_permitted_pipelines_1<aa_pipeline::component_option::zero_dm, false> runner(strategy, input_data.data());
  if(runner.setup()) {
    while(runner.next(true)) {
      std::cout << "NOTICE: Pipeline running over next chunk." << std::endl;
    }
  }
  
  std::cout << "NOTICE: Finished." << std::endl;
  return 0;
}
