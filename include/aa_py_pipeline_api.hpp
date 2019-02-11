#ifndef ASTRO_ACCELERATE_AA_PY_PIPELINE_API_HPP
#define ASTRO_ACCELERATE_AA_PY_PIPELINE_API_HPP

#include "aa_pipeline_api.hpp"
#include "aa_py_filterbank_metadata.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {
      aa_pipeline_api<unsigned short>* aa_py_pipeline_api(const aa_py_filterbank_metadata_struct metadata, unsigned short const*const input_data, const int card_number);
      void aa_py_pipeline_api_delete(aa_pipeline_api<unsigned short> const*const obj);
      bool aa_py_pipeline_api_bind_ddtr_plan(aa_pipeline_api<unsigned short> *const obj, aa_ddtr_plan const*const plan);
      bool aa_py_pipeline_api_bind_analysis_plan(aa_pipeline_api<unsigned short> *const obj, aa_analysis_plan const*const plan)
      bool aa_py_pipeline_api_run(aa_pipeline_api<unsigned short> *const obj);
    }
  }
}

#endif
