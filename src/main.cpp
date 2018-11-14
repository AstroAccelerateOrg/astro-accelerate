//
//  main.cpp
//  aapipeline
//
//  Created by Cees Carels on Monday 22/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_config.hpp"

#include "aa_sigproc_input.hpp"
#include "aa_ddtr_pipeline.hpp"

using namespace astroaccelerate;

int main(int argc, char *argv[]) {
  // If an input_file should be read, then parameters are collected here.
  aa_CLI cli;
  for(int i = 0; i < argc; i++) {
    cli.input.push_back(argv[i]);
  }
  aa_config configuration(cli);
  configuration.setup();


  
  std::cout << "NOTICE: Finished." << std::endl;
  return 0;
}
