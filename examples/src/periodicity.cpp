/**
 * Example code for linking against astro-accelerate library.
 *
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda/samples/common/inc/ -lcudart dedispersion.cpp -o test
 */

#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_permitted_pipelines_generic.hpp"

#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"

#include "aa_periodicity_plan.hpp"
#include "aa_periodicity_strategy.hpp"

using namespace astroaccelerate;

int main() {
//-------------- Select de-dispersion plan
  aa_ddtr_plan ddtr_plan;
  ddtr_plan.add_dm(0, 370, 0.307, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).
  ddtr_plan.add_dm(370, 740, 0.652, 2, 2);
//--------------<

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
  
  int device = 0;
  aa_device_info selected_device(device);

  //-------------- Configure pipeline. Select components and their options
  aa_pipeline::pipeline pipeline_components;
  pipeline_components.insert(aa_pipeline::component::dedispersion); // pipeline must always contain dedispersion step
  pipeline_components.insert(aa_pipeline::component::analysis); //optional
  pipeline_components.insert(aa_pipeline::component::periodicity); // optional
  //pipeline.insert(aa_pipeline::component::fdas); // optional
  
  aa_pipeline::pipeline_option pipeline_options;
  pipeline_options.insert(aa_pipeline::component_option::zero_dm);
  //--------------<

  bool enable_analysis = true;       // The strategy will be optimised to run just dedispersion
  aa_ddtr_strategy ddtr_strategy(ddtr_plan, metadata, selected_device.free_memory(), enable_analysis, &selected_device);
  
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
  const bool  enable_MSD_outlier_rejection = true;
  const aa_analysis_plan::selectable_candidate_algorithm algo = aa_analysis_plan::selectable_candidate_algorithm::peak_find;
  
  aa_analysis_plan analysis_plan(ddtr_strategy, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, algo, enable_MSD_outlier_rejection);
  aa_analysis_strategy analysis_strategy(analysis_plan, &selected_device);

  if(!(analysis_strategy.ready())) {
    std::cout << "ERROR: analysis_strategy not ready." << std::endl;
    return 0;
  }
  
  const float periodicity_sigma_cutoff = 0.0;
  const float periodicity_sigma_constant = sigma_constant;
  const int   nHarmonics = 3;
  const int   export_powers = 0;
  const bool  candidate_algorithm = false;
  const bool  enable_outlier_rejection = false;
  
  aa_periodicity_plan periodicity_plan(periodicity_sigma_cutoff, periodicity_sigma_constant, nHarmonics, export_powers, candidate_algorithm, enable_outlier_rejection);
  aa_periodicity_strategy periodicity_strategy(periodicity_plan);

  if(!periodicity_strategy.ready()) {
    std::cout << "ERROR: periodicity_strategy not ready." << std::endl;
  }

  //-------------- Create empty strategy object for unused components
  aa_fdas_strategy empty_fdas_strategy;
  aa_jerk_strategy empty_jerk_strategy;
  //--------------<

  aa_permitted_pipelines_generic pipeline_runner(pipeline_components, pipeline_options, ddtr_strategy, analysis_strategy, periodicity_strategy, empty_fdas_strategy, empty_jerk_strategy, false, false, false, false, false, input_data.data());  
//  aa_permitted_pipelines_3<aa_pipeline::component_option::zero_dm, false> runner(ddtr_strategy, analysis_strategy, periodicity_strategy, input_data.data());
  if(pipeline_runner.setup()) {
    while(pipeline_runner.next()) {
      std::cout << "NOTICE: Pipeline running over next chunk." << std::endl;
    }
  }
  
  std::cout << "NOTICE: Finished." << std::endl;
  return 0;
}
