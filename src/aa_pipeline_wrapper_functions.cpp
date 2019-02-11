#include "aa_pipeline_wrapper_functions.hpp"

namespace astroaccelerate {

  /** \brief Function that only performs dedispersion, and uses a raw pointer to input data.*/
  void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_pipeline::pipeline_option &pipeline_options, std::vector<aa_ddtr_plan::dm> dm_ranges, unsigned short const*const input_data) {
    std::vector<aa_pipeline::component> selected_components = {aa_pipeline::component::dedispersion};
    aa_pipeline_generic(selected_components, pipeline_options, filterbank_data, dm_ranges, input_data);
  }

  /** \brief Function that only performs dedispersion, and uses a std::vector to input data. */
  void dedisperse_telescope_data(const aa_filterbank_metadata &filterbank_data, const aa_pipeline::pipeline_option &pipeline_options, std::vector<aa_ddtr_plan::dm> dm_ranges, const std::vector<unsigned short> &input_data) {
    // This overloaded function uses a std::vector<unsigned short> and calls the function version using a raw pointer to the data.
    dedisperse_telescope_data(filterbank_data, pipeline_options, dm_ranges, input_data.data());
  }

  /** \brief Function that performs dedispersion and analysis, and uses a raw pointer to input data. */
  void dedisperse_analyse_telescope_data(const aa_filterbank_metadata &filterbank_data,
					 const aa_pipeline::pipeline_option &pipeline_options,
					 std::vector<aa_ddtr_plan::dm> dm_ranges,
					 unsigned short const*const input_data,
					 const float &sigma_cutoff,
					 const float &sigma_constant,
					 const float &max_boxcar_width_in_sec,
					 const bool &enable_candidate_algorithm,
					 const bool &enable_msd_baseline_noise_algorithm) {
    std::vector<aa_pipeline::component> selected_components = {aa_pipeline::component::dedispersion, aa_pipeline::component::analysis};
    aa_pipeline_generic(selected_components, pipeline_options, filterbank_data, dm_ranges, input_data, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, enable_candidate_algorithm, enable_msd_baseline_noise_algorithm);
  }

  /** \brief Function that performs dedispersion and analysis, and uses a std::vector to input data. */
  void dedisperse_analyse_telescope_data(const aa_filterbank_metadata &filterbank_data,
					 const aa_pipeline::pipeline_option &pipeline_options,
					 std::vector<aa_ddtr_plan::dm> dm_ranges,
					 const std::vector<unsigned short> &input_data,
					 const float &sigma_cutoff,
					 const float &sigma_constant,
					 const float &max_boxcar_width_in_sec,
					 const bool &enable_candidate_algorithm,
					 const bool &enable_msd_baseline_noise_algorithm) {
    // This overloaded function uses a std::vector<unsigned short> and calls the function version using a raw pointer to the data.
    dedisperse_analyse_telescope_data(filterbank_data, pipeline_options, dm_ranges, input_data.data(), sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, enable_candidate_algorithm, enable_msd_baseline_noise_algorithm);
  }

} //namespace astroaccelerate
