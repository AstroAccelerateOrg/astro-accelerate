#ifndef ASTRO_ACCELERATE_AA_PIPELINE_GENERIC_HPP
#define ASTRO_ACCELERATE_AA_PIPELINE_GENERIC_HPP

#include <iostream>

#include "aa_sigproc_input.hpp"
#include "aa_pipeline.hpp"
#include "aa_config.hpp"
#include "aa_pipeline_api.hpp"

namespace astroaccelerate {

  /**
   * \brief Templated function that takes a pipeline and pipeline details, and uses the API to process the corresponding pipeline.
   * \details No pipeline will run if the input parameters are invalid.
   * \details This function serves as boilerplate code that provides a wrapper around the API.
   * \details This function serves as an example code for library users to integrate AstroAccelerate into their own applications.
   * \author AstroAccelerate.
   * \date 16 August 2019.
   */  
  template <typename T>
  void aa_pipeline_generic(const std::vector<aa_pipeline::component> &selected_components,
			   const aa_pipeline::pipeline_option &pipeline_options,
			   const aa_filterbank_metadata &filterbank_data,
			   std::vector<aa_ddtr_plan::dm> dm_ranges,
			   T const*const input_data,
			   const float &analysis_sigma_cutoff = 0.0,
			   const float &analysis_sigma_constant = 0.0,
			   const float &analysis_max_boxcar_width_in_sec = 0.0,
			   const bool  &analysis_enable_msd_baseline_noise_algorithm = false,
			   const float &periodicity_sigma_cutoff = 0.0,
			   const float &periodicity_sigma_constant = 0.0,
			   const int   &periodicity_nHarmonics = 0.0,
			   const int   &periodicity_export_powers = 0,
			   const bool  &periodicity_candidate_algorithm = false,
			   const bool  &periodicity_enable_outlier_rejection = false) {
    /**
     * Boilerplate code for executing a pipeline of components.
     * 
     * Default parameters for aa_pipeline::analysis are set.
     * Although this function will configure an aa_analysis_plan in all cases,
     * if std::vector<aa_pipeline::component> &selected_components does not contain
     * aa_pipeline::component::analysis, then aa_pipeline_api will ignore aa_analysis_plan
     * when this function attempts to bind aa_analysis_plan.
     */
    
    // Configure astro-accelerate as a library user
    aa_pipeline::pipeline the_pipeline;
    //the_pipeline = aa_permitted_pipelines::pipeline1;   // EITHER: Use a pre-configured pipeline
    
    //OR insert component manually
    for(size_t i = 0; i < selected_components.size(); i++) {
      the_pipeline.insert(selected_components.at(i));
    }
    
	
	int device = 0;
	aa_device_info selected_device(device);
        
    aa_config configuration(the_pipeline);   // Set the pipeline and other run settings that would come from an input_file
    the_pipeline = configuration.setup();    // The configuration validates whether the pipeline is valid and returns either a valid pipeline or a trivial pipeline
    
    // Supply the requested pipeline and telescope data to a pipeline manager, which will check which components are required to be configured.
    // If a component is not required, then even if it is supplied, it will be ignored.
    aa_pipeline_api<T> pipeline_manager(the_pipeline, pipeline_options, filterbank_data, input_data, selected_device);
    
    // Bind the Plan to the manager
    aa_ddtr_plan ddtr_plan;
    ddtr_plan.add_dm(dm_ranges);

    if(pipeline_manager.bind(ddtr_plan)) {
      LOG(log_level::notice, "ddtr_plan bound successfully.");
    }
    else {
      LOG(log_level::error, "Could not bind ddtr_plan.");
    }

    // When aa_ddtr_plan was bound to aa_pipeline_api, aa_pipeline_api already
    // knew whether aa_pipeline::component::analysis was supplied.
    // Therefore, aa_ddtr_strategy was computed depending on
    // whether analysis would be requested.
    // As such, the aa_ddtr_strategy that will be passed to this method
    // will be correctly configured.
    //
    // Also, since it is not possible to obtain an aa_ddtr_strategy
    // without an aa_ddtr_plan, the user is unable to attempt to configure
    // aa_analysis_plan without first configuring
    // aa_ddtr_plan and obtaining an aa_ddtr_strategy.
    //
    // Therefore, if the user wishes not to use the API, they may manually
    // interact with the interface of the aa_ddtr_plan class and the
    // aa_ddtr_strategy class.
    // In this case, when the user attempts to obtain an aa_ddtr_strategy
    // from an aa_ddtr_plan, they must again specify whether they will
    // require to run analysis.
    // Lastly, aa_ddtr_strategy contains a member field to query whether
    // analysis will be run. This enables aa_analysis_strategy to validate
    // the aa_ddtr_strategy that was supplied to it.
	
	
	aa_analysis_plan::selectable_candidate_algorithm selected_candidate_algorithm;
	if (pipeline_options.find(aa_pipeline::component_option::candidate_algorithm) != pipeline_options.end()) {
		selected_candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::threshold;
	}
	else if (pipeline_options.find(aa_pipeline::component_option::candidate_filtering) != pipeline_options.end()) {
		selected_candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::peak_filtering;
	}
    else {
      selected_candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::peak_find ;
    }
	
    
    aa_analysis_plan analysis_plan(pipeline_manager.ddtr_strategy(),
				   analysis_sigma_cutoff,
				   analysis_sigma_constant,
				   analysis_max_boxcar_width_in_sec,
				   selected_candidate_algorithm,
				   analysis_enable_msd_baseline_noise_algorithm);
    pipeline_manager.bind(analysis_plan);

    aa_periodicity_plan periodicity_plan(periodicity_sigma_cutoff,
					 periodicity_sigma_constant,
					 periodicity_nHarmonics,
					 periodicity_export_powers,
					 periodicity_candidate_algorithm,
					 periodicity_enable_outlier_rejection);
    pipeline_manager.bind(periodicity_plan);
    
    // Bind further plans as necessary
    // ...
    
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
      LOG(log_level::notice, "The pipeline could not start or had errors.");
    }
  }  
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PIPELINE_GENERIC_HPP
