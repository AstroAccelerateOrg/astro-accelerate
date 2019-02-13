#include "aa_py_ddtr_strategy.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {
      void aa_py_ddtr_strategy_delete(aa_ddtr_strategy const*const obj) {
	delete obj;
      }
    } // extern "C"
  } // namespace python
} // namespace astroaccelerate
  
