//
//  aa_ddtr_pipeline.cpp
//  aapipeline
//
//  Created by Cees Carels on Wednesday 24/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_pipeline_wrapper_functions.hpp"

namespace astroaccelerate {

  void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_compute::pipeline_detail &pipeline_details, std::vector<aa_ddtr_plan::dm> dm_ranges, unsigned short const*const input_data, float *&output_data) {
    std::vector<aa_compute::modules> selected_modules = {aa_compute::modules::dedispersion};
    aa_pipeline_generic(selected_modules, pipeline_details, filterbank_data, dm_ranges, input_data, output_data);
  }

  void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_compute::pipeline_detail &pipeline_details, std::vector<aa_ddtr_plan::dm> dm_ranges, const std::vector<unsigned short> &input_data, float *&output_data) {
    // This overloaded function uses a std::vector<unsigned short> and calls the function version using a raw pointer to the data.
    dedisperse_telescope_data(filterbank_data, pipeline_details, dm_ranges, input_data.data(), output_data);
  }

  void dedisperse_analyse_telescope_data(const aa_filterbank_metadata &filterbank_data,
					 const aa_compute::pipeline_detail &pipeline_details,
					 std::vector<aa_ddtr_plan::dm> dm_ranges,
					 unsigned short const*const input_data,
					 float *&output_data,
					 const float &sigma_cutoff,
					 const float &sigma_constant,
					 const float &max_boxcar_width_in_sec,
					 const bool &enable_candidate_algorithm,
					 const bool &enable_sps_baseline_noise_algorithm) {
    std::vector<aa_compute::modules> selected_modules = {aa_compute::modules::dedispersion, aa_compute::modules::analysis};
    aa_pipeline_generic(selected_modules, pipeline_details, filterbank_data, dm_ranges, input_data, output_data, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, enable_candidate_algorithm, enable_sps_baseline_noise_algorithm);
  }
  
  void dedisperse_analyse_telescope_data(const aa_filterbank_metadata &filterbank_data,
					 const aa_compute::pipeline_detail &pipeline_details,
					 std::vector<aa_ddtr_plan::dm> dm_ranges,
					 const std::vector<unsigned short> &input_data,
					 float *&output_data,
					 const float &sigma_cutoff,
					 const float &sigma_constant,
					 const float &max_boxcar_width_in_sec,
					 const bool &enable_candidate_algorithm,
					 const bool &enable_sps_baseline_noise_algorithm) {
    // This overloaded function uses a std::vector<unsigned short> and calls the function version using a raw pointer to the data.
    dedisperse_analyse_telescope_data(filterbank_data, pipeline_details, dm_ranges, input_data.data(), output_data, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, enable_candidate_algorithm, enable_sps_baseline_noise_algorithm);
  }

} //namespace astroaccelerate
