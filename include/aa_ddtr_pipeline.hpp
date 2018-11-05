//
//  aa_ddtr_pipeline.hpp
//  aapipeline
//
//  Created by Cees Carels on Wednesday 24/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_ddtr_pipeline_hpp
#define aa_ddtr_pipeline_hpp

#include <stdio.h>
#include <vector>
#include "aa_compute.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_pipeline_generic.hpp"

void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, unsigned short *input_data, float *output_data);
void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, std::vector<unsigned short> input_data, float *output_data);

#endif /* aa_ddtr_pipeline_hpp */
