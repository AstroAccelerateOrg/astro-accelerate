//
//  aa_ddtr_pipeline.cpp
//  aapipeline
//
//  Created by Cees Carels on Wednesday 24/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_ddtr_pipeline.hpp"

void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, std::vector<aa_ddtr_plan::dm> dm_ranges, unsigned short *input_data, float *output_data) {
    std::vector<aa_compute::modules> selected_modules = {aa_compute::modules::dedispersion};
    aa_pipeline_generic(selected_modules, filterbank_data, dm_ranges, input_data, output_data);
}

void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, std::vector<aa_ddtr_plan::dm> dm_ranges, std::vector<unsigned short> input_data, float *output_data) {
    std::vector<aa_compute::modules> selected_modules = {aa_compute::modules::dedispersion};
    aa_pipeline_generic(selected_modules, filterbank_data, dm_ranges, input_data.data(), output_data);
}
