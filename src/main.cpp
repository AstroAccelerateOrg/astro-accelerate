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
    
      aa_sigproc_input       filterbank_datafile("/mnt/data/AstroAccelerate/filterbank/ska-mid-b2.fil");
      aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();
      filterbank_datafile.read_telescope();
      float *my_output_data = NULL;

      std::vector<aa_ddtr_plan::dm> dm_ranges;
      aa_ddtr_plan::dm tmp1 = {0, 150, 0.1, 1, 1};
      aa_ddtr_plan::dm tmp2 = {150, 300, 0.2, 1, 1};
      dm_ranges.push_back(tmp1);
      dm_ranges.push_back(tmp2);
      
      dedisperse_telescope_data(filterbank_metadata, dm_ranges, filterbank_datafile.input_buffer(), my_output_data);
    
      return 0;
}
