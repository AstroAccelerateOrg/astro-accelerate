//
//  aa_ddtr_pipeline.hpp
//  aapipeline
//
//  Created by Cees Carels on Wednesday 24/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_DDTR_PIPELINE_HPP
#define ASTRO_ACCELERATE_DDTR_PIPELINE_HPP

#include <stdio.h>
#include <vector>
#include "aa_compute.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_pipeline_generic.hpp"

void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, std::vector<aa_ddtr_plan::dm> dm_ranges, unsigned short *input_data, float *&output_data);
void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, std::vector<aa_ddtr_plan::dm> dm_ranges, std::vector<unsigned short> input_data, float *&output_data);

#endif /* ASTRO_ACCELERATE_DDTR_PIPELINE_HPP */
