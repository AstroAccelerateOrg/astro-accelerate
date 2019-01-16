/*! \mainpage AstroAccelerate
 *
 * \section intro_sec Introduction
 *
 * AstroAccelerate is a many core accelerated software package for processing time domain radio-astronomy data.
 *
 * \section install_sec Usage Guide
 * Please refer to the README.md included in the repository for instructions on how to use the standalone version of AstroAccelerate.
 *
 * Please refer to the MANUAL.md included in the repoistory for instructions on how to use the AstroAccelerate API in your own code.
 * 
 * Please refer to the Developers Guide in CONTRIBUTING.md for instructions on how to develop and contribute to the repository.
 *
 * \section contact_sec Contact
 * Please see contact details in the repository and in COLLABORATORS.md.
 * 
 */

#include "aa_config.hpp"
#include "aa_pipeline_api.hpp"
#include "aa_pipeline.hpp"
#include "aa_sigproc_input.hpp"
#include "aa_pipeline_wrapper_functions.hpp"
#include "aa_params.hpp"
#include "aa_device_info.hpp"

#include "aa_log.hpp"

using namespace astroaccelerate;

int main(int argc, char *argv[]) {
  aa_command_line_arguments cli;
  for(int i = 0; i < argc; i++) {
    cli.input.push_back(argv[i]);
  }
  aa_config cli_configuration(cli);

  aa_ddtr_plan ddtr_plan;
  std::string file_path;
  aa_config_flags user_flags = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, std::vector<aa_pipeline::debug>()};
  aa_pipeline::pipeline_option pipeline_options;
  aa_pipeline::pipeline pipeline = cli_configuration.setup(ddtr_plan, user_flags, pipeline_options, file_path);

  LOG(log_level::notice, "File path "+file_path);
  for(auto const i : pipeline) {
    LOG(log_level::notice, component_name(i));
  }
  
  aa_sigproc_input       filterbank_datafile(file_path.c_str());
  aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();

  if(!filterbank_datafile.read_signal()) {
    LOG(log_level::error, "Could not read telescope data.");
    return 0;
  }
  
  //Select card
  aa_device_info& device_info = aa_device_info::instance();
  if(device_info.check_for_devices()) {
    LOG(log_level::notice, "Checked for devices.");
  }
  else {
    LOG(log_level::error, "Could not find any devices.");
  }
  
  aa_device_info::CARD_ID selected_card = CARD;
  aa_device_info::aa_card_info selected_card_info;
  if(device_info.init_card(selected_card, selected_card_info)) {
    LOG(log_level::notice, "init_card complete.");
  }
  else {
    LOG(log_level::error, "init_card incomplete.")
  }
  
  aa_config configuration(pipeline);   // Set the pipeline and other run settings that would come from an input_file
  
  aa_pipeline_api<unsigned short> pipeline_manager(pipeline,
						   pipeline_options,
						   filterbank_metadata,
						   filterbank_datafile.input_buffer().data(),
						   selected_card_info);
  
  if(pipeline_manager.bind(ddtr_plan)) {
    LOG(log_level::notice, "ddtr_plan bound successfully.");
  }
  else {
    LOG(log_level::error, "Could not bind ddtr_plan.");
    //    return 0;
  }
  
  aa_analysis_plan::selectable_candidate_algorithm candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::off;
  if(user_flags.candidate_algorithm) {
    candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::on;
  }

  bool sps_baseline_noise = false;
  if(pipeline_options.find(aa_pipeline::component_option::sps_baseline_noise) != pipeline_options.end()) {
    sps_baseline_noise = true;
  }

  if(pipeline.find(aa_pipeline::component::analysis) != pipeline.end()) {
    aa_analysis_plan analysis_plan(pipeline_manager.ddtr_strategy(),
				   user_flags.sigma_cutoff,
				   user_flags.sigma_constant,
				   user_flags.max_boxcar_width_in_sec,
				   candidate_algorithm,
				   sps_baseline_noise);
    if(pipeline_manager.bind(analysis_plan)) {
      LOG(log_level::notice, "analysis_plan bound successfully.");
    }
    else {
      LOG(log_level::error, "Could not bind analysis_plan.");
    }
  }
  
  if(pipeline.find(aa_pipeline::component::periodicity) != pipeline.end()) {
    //If these settings come from the input_file, then move them into aa_config to be read from the file.
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

  if(pipeline.find(aa_pipeline::component::fdas) != pipeline.end()) {
    aa_fdas_plan fdas_plan(user_flags.sigma_cutoff,
			   user_flags.nboots,
			   user_flags.ntrial_bins,
			   user_flags.navdms,
			   user_flags.narrow,
			   user_flags.wide,
			   user_flags.nsearch,
			   user_flags.aggression);
    pipeline_manager.bind(fdas_plan);
  }
  
  for(size_t i = 0; i < ddtr_plan.range(); i++) {
    LOG(log_level::dev_debug, std::to_string(ddtr_plan.user_dm(i).low)
	+ " " + std::to_string(ddtr_plan.user_dm(i).high)
	+ " " + std::to_string(ddtr_plan.user_dm(i).step)
	+ " " + std::to_string(ddtr_plan.user_dm(i).inBin)
	+ " " + std::to_string(ddtr_plan.user_dm(i).outBin));
  }

  // Validate if all Plans and Strategies are valid and ready to run
  // Optional: Add throw catch to force user to check their settings
  if(pipeline_manager.ready()) {
    LOG(log_level::notice, "Pipeline is ready.");
  }
  else {
    LOG(log_level::notice, "Pipeline is not ready.");
  }
  
  // Run the pipeline
  if(pipeline_manager.run()) {
    LOG(log_level::notice, "The pipeline finished successfully.");
  }
  else {
    LOG(log_level::error, "The pipeline could not start or had errors.");
  }
  
  LOG(log_level::notice, "Finished.");
  return 0;
}
