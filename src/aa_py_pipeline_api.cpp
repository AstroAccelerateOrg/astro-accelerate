#include "aa_py_pipeline_api.hpp"
#include "aa_log.hpp"
#include "aa_permitted_pipelines.hpp"
#include <iostream>

namespace astroaccelerate {
  namespace python {
    extern "C" {
      aa_pipeline_api<unsigned short>* aa_py_pipeline_api(const aa_py_filterbank_metadata_struct aa_py_metadata, unsigned short const*const input_data, const int card_number) {
	LOG(log_level::notice, "Creating aa_py_pipeline_api in library.");
	
	aa_filterbank_metadata metadata(aa_py_metadata.m_tstart, aa_py_metadata.m_tsamp, aa_py_metadata.m_nbits, aa_py_metadata.m_nsamples, aa_py_metadata.m_fch1, aa_py_metadata.m_foff, aa_py_metadata.m_nchans);
	
	aa_device_info& device_info = aa_device_info::instance();
	if(device_info.check_for_devices()) {
	  LOG(log_level::notice, "Checked for devices.");
	}
	else {
	  LOG(log_level::error, "Could not find any devices.");
	}
	
	aa_device_info::CARD_ID selected_card = card_number;
	aa_device_info::aa_card_info selected_card_info;
	if(device_info.init_card(selected_card, selected_card_info)) {
	  LOG(log_level::notice, "init_card complete. Selected card " + std::to_string(selected_card) + ".");
	}
	else {
	  LOG(log_level::error, "init_card incomplete.");
	}
	
	aa_pipeline::pipeline requested_pipeline = aa_permitted_pipelines::pipeline5;
	aa_pipeline::pipeline_option pipeline_options;
	
	aa_pipeline_api<unsigned short>* ptr = new aa_pipeline_api<unsigned short>(requested_pipeline, pipeline_options, metadata, input_data, selected_card_info);
	return ptr;
      }
      
      void aa_py_pipeline_api_delete(aa_pipeline_api<unsigned short> const*const obj) {
	delete obj;
      }

      bool aa_py_pipeline_api_bind_ddtr_plan(aa_pipeline_api<unsigned short> *const obj, aa_ddtr_plan const*const plan) {
	return obj->bind(*plan);
      }

      aa_ddtr_strategy* aa_py_pipeline_api_ddtr_strategy(aa_pipeline_api<unsigned short> *const obj) {
	return new aa_ddtr_strategy(obj->ddtr_strategy());
      }
      
      bool aa_py_pipeline_api_bind_analysis_plan(aa_pipeline_api<unsigned short> *const obj, aa_analysis_plan const*const plan) {
	return obj->bind(*plan);
      }

      aa_analysis_strategy* aa_py_pipeline_api_analysis_strategy(aa_pipeline_api<unsigned short> *const obj) {
	return new aa_analysis_strategy(obj->analysis_strategy());
      }

      bool aa_py_pipeline_api_bind_periodicity_plan(aa_pipeline_api<unsigned short> *const obj, const float sigma_cutoff, const float sigma_constant, const int nHarmonics, const int export_powers, const bool candidate_algorithm, const bool enable_msd_baseline_noise) {
	aa_periodicity_plan plan(sigma_cutoff, sigma_constant, nHarmonics, export_powers, candidate_algorithm, enable_msd_baseline_noise);
	return obj->bind(plan);
      }

      bool aa_py_pipeline_api_bind_fdas_plan(aa_pipeline_api<unsigned short> *const obj, const float sigma_cutoff, const float sigma_constant, const int num_boots, const int num_trial_bins, const int navdms, const float narrow, const float wide, const int nsearch, const float aggression, const bool enable_msd_baseline_noise) {
	aa_fdas_plan plan(sigma_cutoff, sigma_constant, num_boots, num_trial_bins, navdms, narrow, wide, nsearch, aggression, enable_msd_baseline_noise);
	return obj->bind(plan);
      }

      aa_fdas_strategy aa_py_pipeline_api_fdas_strategy(aa_pipeline_api<unsigned short> *const obj) {
	return obj->fdas_strategy();
      }
      
      bool aa_py_pipeline_api_run(aa_pipeline_api<unsigned short> *const obj) {
	std::cout << "run at " << obj << std::endl;
	if(obj->ready()) {
	  return obj->run();
	}
	else {
	  return false;
	}
      }
    }
  }
}
