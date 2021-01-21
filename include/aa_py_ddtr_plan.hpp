#ifndef ASTRO_ACCELERATE_AA_PY_DDTR_PLAN_HPP
#define ASTRO_ACCELERATE_AA_PY_DDTR_PLAN_HPP

#include "aa_ddtr_plan.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {
      aa_ddtr_plan* aa_py_ddtr_plan();
      void aa_py_ddtr_plan_delete(aa_ddtr_plan const*const obj);

      bool aa_py_ddtr_plan_add_dm(aa_ddtr_plan *const obj,
		  const float low,
		  const float high,
		  const float step,
		  const int inBin,
		  const int outBin);
      size_t aa_py_ddtr_plan_range(aa_ddtr_plan const*const obj);

      bool aa_py_ddtr_plan_set_power(aa_ddtr_plan *const obj, const float power);
      float aa_py_ddtr_plan_power(aa_ddtr_plan const*const obj);
      
      bool aa_py_ddtr_plan_set_enable_msd_baseline(aa_ddtr_plan *const obj, const bool flag);
      bool aa_py_ddtr_plan_enable_msd_baseline_noise(aa_ddtr_plan const*const obj);
	  
	  bool aa_py_ddtr_plan_bind_bandpass_normalization(aa_ddtr_plan *const obj, float const*const data, const int data_size);
    }
  }
}

#endif // ASTRO_ACCELERATE_AA_PY_DDTR_PLAN_HPP
