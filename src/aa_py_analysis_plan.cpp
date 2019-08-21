#include "aa_py_analysis_plan.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {
      aa_analysis_plan* aa_py_analysis_plan(aa_ddtr_strategy const*const ddtr_strategy, const float sigma_cutoff, const float sigma_constant, const float max_boxcar_width_in_sec, const int candidate_algorithm, const bool enable_msd_baseline_noise) {
	aa_analysis_plan::selectable_candidate_algorithm candidate_algorithm_flag = aa_analysis_plan::selectable_candidate_algorithm::peak_find;
	if(candidate_algorithm == 1) {
	  candidate_algorithm_flag = aa_analysis_plan::selectable_candidate_algorithm::threshold;
	} else if (candidate_algorithm == 2) {
		candidate_algorithm_flag = aa_analysis_plan::selectable_candidate_algorithm::peak_filtering;
	} else if (candidate_algorithm == 0) {
		candidate_algorithm_flag = aa_analysis_plan::selectable_candidate_algorithm::peak_find;
	}

	return new aa_analysis_plan(*ddtr_strategy, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, candidate_algorithm_flag, enable_msd_baseline_noise);
      }

      void aa_analysis_plan_delete(aa_analysis_plan const*const obj) {
	delete obj;
      }

      float aa_py_analysis_plan_sigma_cutoff(aa_analysis_plan const*const obj) {
	return obj->sigma_cutoff();
      }
      
      float aa_py_analysis_plan_sigma_constant(aa_analysis_plan const*const obj) {
	return obj->sigma_constant();
      }
      
      float aa_py_analysis_plan_max_boxcar_width_in_sec(aa_analysis_plan const*const obj) {
	return obj->max_boxcar_width_in_sec();
      }
      
//      bool aa_py_analysis_plan_candidate_algorithm(aa_analysis_plan const*const obj) {
//	if(obj->candidate_algorithm() == aa_analysis_plan::selectable_candidate_algorithm::threshold) {
//	  return true;
//	}
//	else {
//	  return false;
//	}
//      }
      
      bool aa_py_analysis_plan_enable_msd_baseline_noise(aa_analysis_plan const*const obj) {
	return obj->enable_msd_baseline_noise();
      }
    } // extern "C"
  } // python
} // astroaccelerate
