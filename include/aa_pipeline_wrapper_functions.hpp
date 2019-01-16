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
#include "aa_pipeline.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_pipeline_generic.hpp"

namespace astroaccelerate {

  /**
   * \brief The following functions provide C-style wrapper functions to commonly run pipelines. They all call a boilerplate code that executes the pipeline using the API functionality.
   * \author Cees Carels.
   * \date 24 October 2018.
   */

  /** \brief Function that only performs dedispersion, and uses a raw pointer to input data.*/
  void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_pipeline::pipeline_option &pipeline_options, std::vector<aa_ddtr_plan::dm> dm_ranges, unsigned short const*const input_data);
  
  /** \brief Function that only performs dedispersion, and uses a std::vector to input data. */
  void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_pipeline::pipeline_option &pipeline_options, std::vector<aa_ddtr_plan::dm> dm_ranges, const std::vector<unsigned short> &input_data);

  /** \brief Function that performs dedispersion and analysis, and uses a raw pointer to input data. */
  void dedisperse_analyse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_pipeline::pipeline_option &pipeline_options, std::vector<aa_ddtr_plan::dm> dm_ranges, unsigned short const*const input_data);

  /** \brief Function that performs dedispersion and analysis, and uses a std::vector to input data. */
  void dedisperse_analyse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_pipeline::pipeline_option &pipeline_options, std::vector<aa_ddtr_plan::dm> dm_ranges, const std::vector<unsigned short> &input_data);

  
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PIPELINE_WRAPPER_FUNCTIONS_HPP
