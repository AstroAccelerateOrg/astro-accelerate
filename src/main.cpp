//
//  main.cpp
//  aapipeline
//
//  Created by Cees Carels on Monday 22/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_sigproc_input.hpp"
#include "aa_ddtr_pipeline.hpp"


int main(int argc, const char * argv[]) {
    
      aa_sigproc_input       filterbank_datafile("/home/carels/filterbank/BenMeerKAT.fil");
      aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();
      filterbank_datafile.read_telescope();
      float *my_output_data = NULL;
      dedisperse_telescope_data(filterbank_metadata, filterbank_datafile.input_buffer(), my_output_data);
    
      return 0;
}
