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
#include "aa_host_rfi.hpp"

#include "aa_welcome_notice.hpp"

using namespace astroaccelerate;

int main(int argc, char *argv[]) {
	welcome_notice();
	aa_command_line_arguments cli;
	for (int i = 0; i < argc; i++) {
		cli.input.push_back(argv[i]);
	}
	aa_config cli_configuration(cli);

	aa_ddtr_plan ddtr_plan;
	std::string file_path;
	//aa_config_flags contains values like sigma_cutoff, card_id, but also rfi which should be pipeline option. It also contain vector of user_debug enumerator which should be independent.
	aa_config_flags user_flags = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false, std::vector<aa_pipeline::debug>() };
	//pipeline options is a set of options for the pipeline like: zero_dm, output_dmt, candidate_algorithm which should be handled better
	aa_pipeline::pipeline_option pipeline_options;
	//aa_config takes all argument as reference which is extremely confusing while still returning something. Must change to pointers and return nothing 
	aa_pipeline::pipeline pipeline = cli_configuration.setup(ddtr_plan, user_flags, pipeline_options, file_path);
	//pipeline_options.insert(aa_pipeline::component_option::timelog_export_to_file);
	LOG(log_level::notice, "File path "+file_path);
	LOG(log_level::notice, "-----------------------------------");
	LOG(log_level::notice, "Pipeline components:");
	for (auto const i : pipeline) {
		LOG(log_level::notice, "    " + component_name(i));
	}
	LOG(log_level::notice, "Component options:");
	for (auto const i : pipeline_options) {
		LOG(log_level::notice, "    " + component_option_description(i));
	}
	LOG(log_level::notice, "-----------------------------------");

	aa_sigproc_input       filterbank_datafile(file_path.c_str());
	aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();

	if (!filterbank_datafile.read_signal()) {
		LOG(log_level::error, "Could not read telescope data.");
		return 0;
	}

	//Select card
	int device = 0;
	aa_device_info selected_device(device);

	//Why this is here we have already configured the pipeline? Not used later delete?
	//aa_config configuration(pipeline);   // Set the pipeline and other run settings that would come from an input_file

	// Move this to pipeline
	if (user_flags.rfi == 1) {
		LOG(log_level::notice, "Performing host RFI reduction. This feature is experimental.");
		rfi(filterbank_metadata.nsamples(), filterbank_metadata.nchans(), filterbank_datafile.input_buffer_modifiable());
	}

	// What is this doing
	aa_pipeline_api<unsigned short> pipeline_manager(
		pipeline, // list of components
		pipeline_options,
		filterbank_metadata,
		filterbank_datafile.input_buffer().data(),
		selected_device);

	// MSD baseline noise should be moved to new component which would be candidate selection
	bool msd_baseline_noise = false;
	if (pipeline_options.find(aa_pipeline::component_option::msd_baseline_noise) != pipeline_options.end()) {
		msd_baseline_noise = true;
	}

	ddtr_plan.set_enable_msd_baseline_noise(msd_baseline_noise);

	if (pipeline_manager.bind(ddtr_plan)) {
		LOG(log_level::notice, "ddtr_plan bound successfully.");
	}
	else {
		LOG(log_level::error, "Could not bind ddtr_plan.");
		//    return 0;
	}

	// Same here move this into independent component
	aa_analysis_plan::selectable_candidate_algorithm candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::peak_find;
	if (user_flags.candidate_algorithm == 1) {
		candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::threshold;
	}
	else if (user_flags.candidate_algorithm == 2) {
		candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::peak_filtering;
	}
	
	if (pipeline.find(aa_pipeline::component::analysis) != pipeline.end()) {
		aa_analysis_plan analysis_plan(
			pipeline_manager.ddtr_strategy(),
			user_flags.sigma_cutoff,
			user_flags.sigma_constant,
			user_flags.max_boxcar_width_in_sec,
			candidate_algorithm,
			msd_baseline_noise);
		if (pipeline_manager.bind(analysis_plan)) {
			LOG(log_level::notice, "analysis_plan bound successfully.");
		}
		else {
			LOG(log_level::error, "Could not bind analysis_plan.");
		}
	}

	
	if (pipeline.find(aa_pipeline::component::periodicity) != pipeline.end()) {
		//If these settings come from the input_file, then move them into aa_config to be read from the file.
		aa_periodicity_plan periodicity_plan(
			user_flags.periodicity_sigma_cutoff,
			user_flags.sigma_constant,
			user_flags.periodicity_nHarmonics,
			user_flags.power,
			user_flags.candidate_algorithm,
			msd_baseline_noise);
		pipeline_manager.bind(periodicity_plan);
	}

	
	if (pipeline.find(aa_pipeline::component::fdas) != pipeline.end()) {
		aa_fdas_plan fdas_plan(
			user_flags.sigma_cutoff,
			user_flags.sigma_constant,
			msd_baseline_noise);
		pipeline_manager.bind(fdas_plan);
	}
	
	if (pipeline.find(aa_pipeline::component::jerk) != pipeline.end()) {
		bool high_precision_filters = false;
		bool interbinning = false;
		bool always_choose_next_power_of_2 = false;
		bool spectrum_whitening = false;
		size_t free_mem, total_mem;
		cudaMemGetInfo(&free_mem,&total_mem);
		free_mem = free_mem*0.90;
		aa_jerk_plan jerk_plan(
			pipeline_manager.ddtr_strategy().nProcessedTimesamples(),
			pipeline_manager.ddtr_strategy().max_ndms(),
			free_mem,
			user_flags.z_max,
			user_flags.z_step,
			user_flags.w_max,
			user_flags.w_step,
			interbinning,
			high_precision_filters,
			user_flags.sigma_cutoff,
			user_flags.periodicity_nHarmonics,
			msd_baseline_noise,
			user_flags.sigma_constant,
			always_choose_next_power_of_2,
			spectrum_whitening
		);
		pipeline_manager.bind(jerk_plan);
	}

	for (size_t i = 0; i < ddtr_plan.range(); i++) {
		LOG(log_level::dev_debug, std::to_string(ddtr_plan.user_dm(i).low)
			+ " " + std::to_string(ddtr_plan.user_dm(i).high)
			+ " " + std::to_string(ddtr_plan.user_dm(i).step)
			+ " " + std::to_string(ddtr_plan.user_dm(i).inBin)
			+ " " + std::to_string(ddtr_plan.user_dm(i).outBin));
	}

	// Way to general pipeline is probably to send list of required components aka 'pipeline' 
	
	
	
	// Validate if all Plans and Strategies are valid and ready to run
	// Optional: Add throw catch to force user to check their settings
	if (pipeline_manager.ready()) {
		LOG(log_level::notice, "Pipeline is ready.");
	}
	else {
		LOG(log_level::notice, "Pipeline is not ready.");
	}

	// Run the pipeline
	if (pipeline_manager.run()) {
		LOG(log_level::notice, "The pipeline finished successfully.");
	}
	else {
		LOG(log_level::error, "The pipeline could not start or had errors.");
	}

	LOG(log_level::notice, "Finished.");
	return 0;
}

// TODO: rename 'pipeline' to 'pipeline_components'
