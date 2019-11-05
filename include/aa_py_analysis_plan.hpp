#ifndef ASTRO_ACCELERATE_AA_PY_ANALYSIS_PLAN_HPP
#define ASTRO_ACCELERATE_AA_PY_ANALYSIS_PLAN_HPP

#include "aa_analysis_plan.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {
      aa_analysis_plan* aa_py_analysis_plan(aa_ddtr_strategy const*const ddtr_strategy, const float sigma_cutoff, const float sigma_constant, const float max_boxcar_width_in_sec, const int candidate_algorithm, const bool enable_msd_baseline_noise);
      void aa_analysis_plan_delete(aa_analysis_plan const*const obj);
      float aa_py_analysis_plan_sigma_cutoff(aa_analysis_plan const*const obj);
      float aa_py_analysis_plan_sigma_constant(aa_analysis_plan const*const obj);
      float aa_py_analysis_plan_max_boxcar_width_in_sec(aa_analysis_plan const*const obj);
      bool aa_py_analysis_plan_candidate_algorithm(aa_analysis_plan const*const obj);
      bool aa_py_analysis_plan_enable_msd_baseline_noise(aa_analysis_plan const*const obj); 
    } // extern "C"
  } // python
} // astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PY_ANALYSIS_PLAN_HPP
