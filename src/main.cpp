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
  // If an input_file should be read, then parameters are collected here.
  aa_CLI cli;
  for(int i = 0; i < argc; i++) {
    cli.input.push_back(argv[i]);
  }
  aa_config configuration(cli);

  aa_ddtr_plan ddtr_plan;
  std::string file_path;
  aa_compute::pipeline pipeline = configuration.setup(ddtr_plan, file_path);

  std::cout << "File path " << file_path << std::endl; 
  
  for(size_t i = 0; i < ddtr_plan.range(); i++) {
    std::cout << ddtr_plan.user_dm(i).low << " " << ddtr_plan.user_dm(i).high << " " << ddtr_plan.user_dm(i).step << " " << ddtr_plan.user_dm(i).inBin << " " << ddtr_plan.user_dm(i).outBin << std::endl;
  }


  for(auto const i : pipeline) {
    std::cout << module_name(i) << std::endl;
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
  
  std::cout << "Going to init card " << std::endl;
    if(device_info.init_card(selected_card, selected_card_info)) {
      std::cout << "NOTICE: init_card complete." << std::endl;
    }
    else {
      std::cout << "ERROR: init_card incomplete." << std::endl;
    }

  
  std::cout << "NOTICE: Finished." << std::endl;
  return 0;
}
