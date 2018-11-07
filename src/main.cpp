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
      if(!filterbank_datafile.read_telescope()) {
	std::cout << "ERROR: Could not read telescope data." << std::endl;
	return 0;
      }
      float *my_output_data = NULL;

      std::vector<aa_ddtr_plan::dm> dm_ranges;
      aa_ddtr_plan::dm range1 = {0, 150, 0.1, 1, 1};
      aa_ddtr_plan::dm range2 = {150, 300, 0.2, 1, 1};
      aa_ddtr_plan::dm range3 = {300, 500, 0.25, 1, 1};
      aa_ddtr_plan::dm range4 = {500, 900, 0.4, 2, 2};
      aa_ddtr_plan::dm range5 = {900, 1200, 0.6, 4, 4};
      aa_ddtr_plan::dm range6 = {1200, 1500, 0.8, 4, 4};
      aa_ddtr_plan::dm range7 = {1500, 2000, 1.0, 4, 4};
      aa_ddtr_plan::dm range8 = {2000, 3000, 2.0, 8, 8};

      dm_ranges.push_back(range1);
      dm_ranges.push_back(range2);
      dm_ranges.push_back(range3);
      dm_ranges.push_back(range4);
      dm_ranges.push_back(range5);
      dm_ranges.push_back(range6);
      dm_ranges.push_back(range7);
      dm_ranges.push_back(range8);
      
      dedisperse_telescope_data(filterbank_metadata, dm_ranges, filterbank_datafile.input_buffer(), my_output_data);
    
      return 0;
}
