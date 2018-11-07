//
//  aa_pipeline_generic.hpp
//  aapipeline
//
//  Created by Cees Carels on Wednesday 24/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_PIPELINE_GENERIC_HPP
#define ASTRO_ACCELERATE_PIPELINE_GENERIC_HPP

#include <iostream>

#include "aa_sigproc_input.hpp"
#include "aa_compute.hpp"
#include "aa_config.hpp"
#include "aa_pipeline.hpp"

template <typename T, typename U>
void aa_pipeline_generic(const std::vector<aa_compute::modules> &selected_modules, const aa_filterbank_metadata &filterbank_data, std::vector<aa_ddtr_plan::dm> dm_ranges, T *input_data, U *output_data) {
    /**
     * Boilerplate code for executing a pipeline of modules
     */
    
    // Configure astro-accelerate as a library user
    aa_compute::pipeline the_pipeline;
    //the_pipeline = aa_permitted_pipelines::pipeline1;   // EITHER: Use a pre-configured pipeline
    
    //OR insert modules manually
    for(size_t i = 0; i < selected_modules.size(); i++) {
        the_pipeline.insert(selected_modules.at(i));
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
        
    aa_config configuration(the_pipeline);   // Set the pipeline and other run settings that would come from an input_file
    the_pipeline = configuration.setup();    // The configuration validates whether the pipeline is valid and returns either a valid pipeline or a trivial pipeline
    
    // Supply the requested pipeline and telescope data to a pipeline manager, which will check which modules are required to be configured
    aa_pipeline<T, U> pipeline_manager(the_pipeline, filterbank_data, selected_card_info);
    
    // Bind the Plan to the manager
    aa_ddtr_plan ddtr_plan;
    for(auto i : dm_ranges) {
      ddtr_plan.add_dm(i);
    }

    if(pipeline_manager.bind(ddtr_plan)) {
        std::cout << "NOTICE: ddtr_plan bound successfully." << std::endl;
    }
    else {
        std::cout << "ERROR: Could not bind ddtr_plan." << std::endl;
    }
    
    aa_analysis_plan analysis_plan;
    pipeline_manager.bind(analysis_plan);
    
    aa_periodicity_plan periodicity_plan;
    pipeline_manager.bind(periodicity_plan);
    
    // Bind further plans as necessary
    // ...
    
    // Bind data
    if(pipeline_manager.bind_data(input_data)) {
        std::cout << "NOTICE: The data was bound to the pipeline successfully." << std::endl;
    }
    else {
        std::cout << "ERROR: The data could not be bound to the pipeline." << std::endl;
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
        std::cout << "NOTICE: Pipeline is not ready" << std::endl;
    }
    
    // Run the pipeline
    if(pipeline_manager.run()) {
        std::cout << "NOTICE: The pipeline finished successfully" << std::endl;
    }
    else {
        std::cout << "NOTICE: The pipeline could not start or had errors" << std::endl;
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
}


#endif /* ASTRO_ACCELERATE_PIPELINE_GENERIC_HPP */
