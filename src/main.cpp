//
//  main.cpp
//  aapipeline
//
//  Created by Cees Carels on Monday 22/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include <iostream>

#include "aa_sigproc_input.hpp"
#include "aa_compute.hpp"
#include "aa_config.hpp"

#include "aa_pipeline.hpp"

#include "aa_ddtr_pipeline.hpp"

#include "aa_sigproc_input.hpp"

int main(int argc, const char * argv[]) {
    
      aa_sigproc_input       filterbank_datafile("/home/carels/filterbank/BenMeerKAT.fil");
      aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();
      filterbank_datafile.read_telescope();
      float *my_output_data = NULL;
      dedisperse_telescope_data(filterbank_metadata, filterbank_datafile.input_buffer(), my_output_data);
    
      return 0;
    
    
    
    /*
     // If an input_file should be read, then parameters are collected here.
     aa_CLI cli;
     for(int i = 0; i < argc; i++) {
     cli.input.push_back(argv[i]);
     }
     */
    //Configure all modules here...
    //Configure through Plan -> Strategy idiom, passing them to the aa_config as a "bind"
    //via a std::move
    
    /*
     //Configure astro-accelerate from an input_file
     const std::string config_file_path = "~/Desktop/input_files/BenMeerKAT.txt";
     aa_config config(config_file_path.c_str(), cli);
     */
    
    //An aa_config can have only one aa_filterbank_metadata associated to it at a time.
    //An aa_config can have only one input data file associated to it at a time.
    //The aa_config will return a trivial Strategy if either the filterbank_metadata *and* input file are not associated to aa_config.
    //If the pipeline is known, and both the input data file and filterbank_metadata are known, and if the Plan settings are valid, then the aa_config will return a non-trivial Strategy.
    //Pass Strategy objects back to aa_config to check all configurations correct.
    //Possibly: implement a unique handle on the Strategy that identifies it to the aa_config.
    
    //The compute pipeline will not be run unless all Strategy objects for the pipeline are valid.
    
    //Passing a new filterbank metadata or input data file to an aa_config object will invalidate all
    //associated Plan and Strategy objects
    
    //However, a user may configure a new Plan, obtain a new Strategy, and run the Pipeline again, without invalidating existing configurations.
    
    
    
    /*std::cout << "Pipeline size " << my_pipeline.size() << std::endl;
    std::cout << aa_compute::module_name(*my_pipeline.begin()) << std::endl;
    std::cout << aa_compute::module_name(*(++my_pipeline.begin())) << std::endl;*/
}
