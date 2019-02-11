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

	return new aa_pipeline_api<unsigned short>(requested_pipeline, pipeline_options, metadata, input_data, selected_card_info);
      }
      
      void aa_py_pipeline_api_delete(aa_pipeline_api<unsigned short> const*const obj) {
	delete obj;
      }

      bool aa_py_pipeline_api_bind_ddtr_plan(aa_pipeline_api<unsigned short> *const obj, aa_ddtr_plan const*const plan) {
	obj->bind(*plan);
	return true;
      }

      bool aa_py_pipeline_api_bind_analysis_plan(aa_pipeline_api<unsigned short> *const obj, aa_analysis_plan const*const plan) {
	obj->bind(*plan);
      }

      bool aa_py_pipeline_api_run(aa_pipeline_api<unsigned short> *const obj) {
	return obj->run();
      }
    }
  }
}
