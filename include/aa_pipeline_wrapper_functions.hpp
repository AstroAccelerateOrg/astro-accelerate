//
//  aa_ddtr_pipeline.hpp
//  aapipeline
//
//  Created by Cees Carels on Wednesday 24/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_AA_PIPELINE_WRAPPER_FUNCTIONS_HPP
#define ASTRO_ACCELERATE_AA_PIPELINE_WRAPPER_FUNCTIONS_HPP

#include <stdio.h>
#include <vector>
#include "aa_compute.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_pipeline_generic.hpp"

namespace astroaccelerate {

  void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_compute::pipeline_detail &pipeline_details, std::vector<aa_ddtr_plan::dm> dm_ranges, unsigned short const*const input_data, float *&output_data);
  void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_compute::pipeline_detail &pipeline_details, std::vector<aa_ddtr_plan::dm> dm_ranges, const std::vector<unsigned short> &input_data, float *&output_data);

  void dedisperse_analyse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_compute::pipeline_detail &pipeline_details, std::vector<aa_ddtr_plan::dm> dm_ranges, unsigned short const*const input_data, float *&output_data);
  void dedisperse_analyse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_compute::pipeline_detail &pipeline_details, std::vector<aa_ddtr_plan::dm> dm_ranges, const std::vector<unsigned short> &input_data, float *&output_data);

  
} //namespace astroaccelerate

#endif /* ASTRO_ACCELERATE_AA_PIPELINE_WRAPPER_FUNCTIONS_HPP */
