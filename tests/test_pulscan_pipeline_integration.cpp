#include <cmath>
#include <iostream>
#include <vector>


#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_pipeline_api.hpp"
#include "aa_permitted_pipelines_generic.hpp"
#include "aa_device_info.hpp"

using namespace astroaccelerate;

int main() {
  std::cout << "Running test_pulscan_pipeline_integration.cpp" << std::endl;

  aa_ddtr_plan ddtr_plan;
  const float dm_low = 0.0f;
  const float dm_high = 10.0f;
  // Keep the shared-memory line within bounds so this exercises the fast
  // dedispersion kernel rather than its automatic failsafe fallback.
  const float dm_step = 0.25f;
  const int in_bin = 1;
  const int out_bin = 1;
  ddtr_plan.add_dm(dm_low, dm_high, dm_step, in_bin, out_bin);

  const double tstart = 58000.0;
  const double tsamp = 6.4e-5;
  const unsigned int nsamples = 4096;
  const int nchans = 64;
  const double fch1 = 1400.0;
  const double bandwidth = 200.0;
  const double foff = -bandwidth / nchans;
  const double nbits = 8;

  aa_filterbank_metadata metadata(tstart, tsamp, nbits, nsamples, fch1, foff, nchans);
  const std::size_t total_samples = static_cast<std::size_t>(metadata.nsamples());
  const std::size_t total_channels = static_cast<std::size_t>(metadata.nchans());

  aa_device_info device(0);
  const size_t free_memory = device.free_memory();

  bool enable_analysis = false;
  aa_ddtr_strategy ddtr_strategy(ddtr_plan, metadata, free_memory, enable_analysis, &device);
  if(!ddtr_strategy.ready()) {
    std::cout << "DDTR strategy not ready" << std::endl;
    return 0;
  }

  std::vector<unsigned short> input(total_samples * total_channels, 0);
  for(std::size_t chan = 0; chan < total_channels; ++chan) {
    for(std::size_t sample = 0; sample < total_samples; ++sample) {
      double value =
          100.0 + 5.0 * std::sin(0.01 * static_cast<double>(sample) +
                                 0.1 * static_cast<double>(chan));
      input[chan * total_samples + sample] = static_cast<unsigned short>(
          std::max(0.0, std::min(255.0, value)));
    }
  }

  aa_pipeline::pipeline components;
  components.insert(aa_pipeline::component::dedispersion);
  components.insert(aa_pipeline::component::periodicity);

  aa_pipeline::pipeline_option options;
  options.insert(aa_pipeline::component_option::copy_ddtr_data_to_host);

  aa_pipeline_api<unsigned short> runner(components, options, metadata, input.data(), device);
  runner.bind(ddtr_plan);

  auto ddtr_strategy_instance = runner.ddtr_strategy();
  if(!ddtr_strategy_instance.ready()) {
    std::cout << "DDTR strategy not ready after binding" << std::endl;
    return 0;
  }

  aa_periodicity_plan periodicity_plan(ddtr_strategy_instance,
                                       10.0f,
                                       true,
                                       3.0f,
                                       32,
                                       PSR_HRMS_GREEDY);
  runner.bind(periodicity_plan);

  if(!runner.ready()) {
    std::cout << "Pipeline runner not ready" << std::endl;
    return 0;
  }

  aa_pipeline_runner::status status_code;
  while(runner.run(status_code)) {
    if(status_code == aa_pipeline_runner::status::error) {
      std::cout << "Pipeline reported error status" << std::endl;
      return 1;
    }
  }

  const auto& candidates = runner.get_pulscan_candidates();
  std::cout << "Pulscan candidate count: " << candidates.size() << std::endl;

  for(std::size_t i = 1; i < candidates.size(); ++i) {
    if(candidates[i - 1].sigma < candidates[i].sigma) {
      std::cout << "Candidate sigma ordering violated" << std::endl;
      return 1;
    }
  }

  std::cout << "Runs" << std::endl;
  return 0;
}
