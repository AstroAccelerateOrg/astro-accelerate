#ifndef ASTRO_ACCELERATE_AA_PY_DDTR_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_PY_DDTR_STRATEGY_HPP

#include "aa_ddtr_strategy.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {
      void aa_py_ddtr_strategy_delete(aa_ddtr_strategy const*const obj);
    } // extern "C"
  } // namespace python
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PY_DDTR_STRATEGY_HPP
