//
//  main.cpp
//  aapipeline
//
//  Created by Cees Carels on Monday 22/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_config.hpp"
#include "aa_pipeline.hpp"
#include "aa_compute.hpp"
#include "aa_sigproc_input.hpp"
#include "aa_pipeline_wrapper_functions.hpp"

#include "aa_device_info.hpp"

using namespace astroaccelerate;

int main(int argc, char *argv[]) {
  aa_CLI cli;
  for(int i = 0; i < argc; i++) {
    cli.input.push_back(argv[i]);
  }
  aa_config cli_configuration(cli);

  aa_ddtr_plan ddtr_plan;
  std::string file_path;
  aa_config_flags user_flags = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, std::vector<aa_compute::debug>()};
  aa_compute::pipeline_detail pipeline_details;
  aa_compute::pipeline pipeline = cli_configuration.setup(ddtr_plan, user_flags, pipeline_details, file_path);

  std::cout << "File path " << file_path << std::endl;
  for(auto const i : pipeline) {
    std::cout << module_name(i) << std::endl;
  }
  
  aa_sigproc_input       filterbank_datafile(file_path.c_str());
  aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();

  if(!filterbank_datafile.read_telescope()) {
    std::cout << "ERROR: Could not read telescope data." << std::endl;
    return 0;
  }
  
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
  
  aa_config configuration(pipeline);   // Set the pipeline and other run settings that would come from an input_file
  float *output_data = NULL;
  aa_pipeline<unsigned short, float> pipeline_manager(pipeline,
						      filterbank_metadata,
						      filterbank_datafile.input_buffer().data(),
						      selected_card_info);
  
  if(pipeline_manager.bind(ddtr_plan)) {
    std::cout << "NOTICE: ddtr_plan bound successfully." << std::endl;
  }
  else {
    std::cout << "ERROR: Could not bind ddtr_plan." << std::endl;
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
    aa_analysis_plan analysis_plan(pipeline_manager.ddtr_strategy(),
				   user_flags.sigma_cutoff,
				   user_flags.sigma_constant,
				   user_flags.max_boxcar_width_in_sec,
				   candidate_algorithm,
				   sps_baseline_noise);
    if(pipeline_manager.bind(analysis_plan)) {
      std::cout << "NOTICE: analysis_plan bound successfully." << std::endl;
    }
    else {
      std::cout << "ERROR: Could not bind analysis_plan." << std::endl;
    }
  }
  
  if(pipeline.find(aa_compute::modules::periodicity) != pipeline.end()) { 
    const float OR_sigma_multiplier = 1.0;              //Is this setting in the input_file? Is it the same one as for analysis?
    const bool periodicity_candidate_algorithm = false; //Is this setting in the input_file? Is it the same one as for analysis?
    const bool enable_outlier_rejection = false;        //Is this setting in the input_file? Is it the same one as for analysis?
    aa_periodicity_plan periodicity_plan(user_flags.sigma_cutoff,
					 OR_sigma_multiplier,
					 user_flags.periodicity_nHarmonics,
					 user_flags.power,
					 periodicity_candidate_algorithm,
					 enable_outlier_rejection);
    
    pipeline_manager.bind(periodicity_plan);
  }
  
  for(size_t i = 0; i < ddtr_plan.range(); i++) {
    std::cout << ddtr_plan.user_dm(i).low << " "
	      << ddtr_plan.user_dm(i).high << " "
	      << ddtr_plan.user_dm(i).step << " "
	      << ddtr_plan.user_dm(i).inBin << " "
	      << ddtr_plan.user_dm(i).outBin << std::endl;
  }

  if(pipeline_manager.transfer_data_to_device()) {
    std::cout << "NOTICE: The data was transferred to the device successfully." << std::endl;
  }
  else {
    std::cout << "ERROR: The data could not be transferred to the device." << std::endl;
  }
  
  // Validate if all Plans and Strategies are valid and ready to run
  // Optional: Add throw catch to force user to check their settings
  if(pipeline_manager.ready()) {
    std::cout << "NOTICE: Pipeline is ready." << std::endl;
  }
  else {
    std::cout << "NOTICE: Pipeline is not ready." << std::endl;
  }
  
  // Run the pipeline
  if(pipeline_manager.run()) {
    std::cout << "NOTICE: The pipeline finished successfully." << std::endl;
  }
  else {
    std::cout << "NOTICE: The pipeline could not start or had errors." << std::endl;
  }

  // Bring data back from device to host                                                      
  if(pipeline_manager.transfer_data_to_host(output_data)) {
    std::cout << "NOTICE: Data was transferred back to host successfully." << std::endl;
  }
  else {
    std::cout << "NOTICE: Data was not transferred back to host." << std::endl;
  }
  
  if(pipeline_manager.unbind_data()) {
    std::cout << "NOTICE: Data was unbound successfully." << std::endl;
  }
  else {
    std::cout << "NOTICE: Data could not be unbound." << std::endl;
  }
  
  std::cout << "NOTICE: Finished." << std::endl;
  return 0;
}
