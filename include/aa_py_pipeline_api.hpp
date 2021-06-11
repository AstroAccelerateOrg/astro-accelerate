#ifndef ASTRO_ACCELERATE_AA_PY_PIPELINE_API_HPP
#define ASTRO_ACCELERATE_AA_PY_PIPELINE_API_HPP

#include "aa_pipeline_api.hpp"
#include "aa_py_filterbank_metadata.hpp"

namespace astroaccelerate {
  /**
   * \namespace astroaccelerate::python
   * \brief A project-wide namespace python is used to separate the implementation of the library from the Python wrapper code.
   */
  namespace python {
    extern "C" {

      /**
       * \struct pipeline_components_struct
       * \brief Wrapper for aa_pipeline::component.
       * \details Please see include/aa_pipeline.hpp for implementation details.
       */
      struct pipeline_components_struct {
	bool dedispersion;
	bool analysis;
	bool periodicity;
	bool fdas;
      };

      /**
       * \struct pipeline_component_option_struct
       * \brief Wrapper for aa_pipeline::component_option.
       * \details Please see include/aa_pipeline.hpp for implementation details.
       */
      struct pipeline_component_option_struct {
	bool zero_dm;
	bool zero_dm_with_outliers;
	bool input_DDTR_normalization;
	bool output_DDTR_normalization;
	bool set_bandpass_average;
	bool old_rfi;
	bool copy_ddtr_data_to_host;
	bool msd_baseline_noise;
	bool output_dmt;
	bool output_ffdot_plan;
	bool output_fdas_list;
	bool candidate_algorithm;
	bool fdas_custom_fft;
	bool fdas_inbin;
	bool fdas_norm;
      };
      
      aa_pipeline_api<unsigned short>* aa_py_pipeline_api(const pipeline_components_struct pipeline, const pipeline_component_option_struct options, const aa_py_filterbank_metadata_struct metadata, unsigned short const*const input_data, const int card_number);
      void aa_py_pipeline_api_delete(aa_pipeline_api<unsigned short> const*const obj);
      bool aa_py_pipeline_api_bind_ddtr_plan(aa_pipeline_api<unsigned short> *const obj, aa_ddtr_plan const*const plan);

      aa_ddtr_strategy* aa_py_pipeline_api_ddtr_strategy(aa_pipeline_api<unsigned short> *const obj);
      
      bool aa_py_pipeline_api_bind_analysis_plan(aa_pipeline_api<unsigned short> *const obj, aa_analysis_plan const*const plan);
      aa_analysis_strategy* aa_py_pipeline_api_analysis_strategy(aa_pipeline_api<unsigned short> *const obj);

      bool aa_py_pipeline_api_bind_periodicity_plan(aa_pipeline_api<unsigned short> *const obj, const float sigma_cutoff, const float sigma_constant, const int nHarmonics, const int export_powers, const bool candidate_algorithm, const bool enable_msd_baseline_noise);

      bool aa_py_pipeline_api_bind_fdas_plan(aa_pipeline_api<unsigned short> *const obj, const float sigma_cutoff, const float sigma_constant, const bool enable_msd_baseline_noise);

      aa_fdas_strategy aa_py_pipeline_api_fdas_strategy(aa_pipeline_api<unsigned short> *const obj);
      
      bool aa_py_pipeline_api_run(aa_pipeline_api<unsigned short> *const obj, int &status_code);
    }
  }
}

#endif
