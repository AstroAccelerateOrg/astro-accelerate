//
//  main.cpp
//  aapipeline
//
//  Created by Cees Carels on Monday 22/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_config.hpp"

#include "aa_sigproc_input.hpp"
#include "aa_pipeline_wrapper_functions.hpp"

#include "aa_device_info.hpp"

using namespace astroaccelerate;

int main(int argc, char *argv[]) {
  aa_CLI cli;
  for(int i = 0; i < argc; i++) {
    cli.input.push_back(argv[i]);
  }
  aa_config configuration(cli);

  aa_ddtr_plan ddtr_plan;
  std::string file_path;
  aa_config_flags user_flags = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, std::vector<aa_compute::debug>()};
  aa_compute::pipeline_detail pipeline_details;
  aa_compute::pipeline pipeline = configuration.setup(ddtr_plan, user_flags, pipeline_details, file_path);

  std::cout << "File path " << file_path << std::endl;
  for(auto const i : pipeline) {
    std::cout << module_name(i) << std::endl;
  }
  
  aa_sigproc_input       filterbank_datafile(file_path.c_str());
  aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();
  
  //Select card
  aa_device_info device_info;
  if(device_info.check_for_devices()) {
    std::cout << "NOTICE: Checked for devices." << std::endl;
  }
  else {
    std::cout << "ERROR: Could not find any devices." << std::endl;
  }
  
  aa_device_info::CARD_ID selected_card = 0;
  aa_device_info::aa_card_info selected_card_info;
  if(device_info.init_card(selected_card, selected_card_info)) {
    std::cout << "NOTICE: init_card complete." << std::endl;
  }
  else {
    std::cout << "ERROR: init_card incomplete." << std::endl;
  }
  
  const size_t free_memory = selected_card_info.free_memory;
  
  bool enable_analysis = false;
  if(pipeline.find(aa_compute::modules::analysis) != pipeline.end()) {
    enable_analysis = true;
  }
    
  aa_ddtr_strategy ddtr_strategy(ddtr_plan, filterbank_metadata, free_memory, enable_analysis);
  if(!ddtr_strategy.ready()) {
    std::cout << "ERROR: Could not configure ddtr_strategy." << std::endl;
    return 0;
  }

  aa_analysis_plan::selectable_candidate_algorithm candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::off;
  if(user_flags.candidate_algorithm) {
    candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::on;
  }

  bool sps_baseline_noise = false;
  if(pipeline_details.find(aa_compute::module_option::sps_baseline_noise) != pipeline_details.end()) {
    sps_baseline_noise = true;
  }

  if(pipeline.find(aa_compute::modules::analysis) != pipeline.end()) {
    aa_analysis_plan analysis_plan(ddtr_strategy,
				   user_flags.sigma_cutoff,
				   user_flags.sigma_constant,
				   user_flags.max_boxcar_width_in_sec,
				   candidate_algorithm,
				   sps_baseline_noise);
    aa_analysis_strategy analysis_strategy(analysis_plan);
    if(!analysis_strategy.ready()) {
      std::cout << "Could not configure analysis_strategy." << std::endl;
      return 0;
    }
  }
  
  if(pipeline.find(aa_compute::modules::periodicity) != pipeline.end()) { 
    const float OR_sigma_multiplier = 1.0;
    const bool periodicity_candidate_algorithm = false;
    const bool enable_outlier_rejection = false;
    aa_periodicity_plan periodicity_plan(user_flags.sigma_cutoff,
					 OR_sigma_multiplier,
					 user_flags.periodicity_nHarmonics,
					 user_flags.power,
					 periodicity_candidate_algorithm,
					 enable_outlier_rejection);
    
    aa_periodicity_strategy periodicity_strategy(periodicity_plan);
    if(!periodicity_strategy.ready()) {
      std::cout << "ERROR: Could not configure periodicity_strategy." << std::endl;
      return 0;
    }
  }
  
  for(size_t i = 0; i < ddtr_plan.range(); i++) {
    std::cout << ddtr_plan.user_dm(i).low << " "
	      << ddtr_plan.user_dm(i).high << " "
	      << ddtr_plan.user_dm(i).step << " "
	      << ddtr_plan.user_dm(i).inBin << " "
	      << ddtr_plan.user_dm(i).outBin << std::endl;
  }

  std::cout << "NOTICE: Finished." << std::endl;
  return 0;
}
