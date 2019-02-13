#ifndef ASTRO_ACCELERATE_AA_PY_SIGPROC_INPUT_HPP
#define ASTRO_ACCELERATE_AA_PY_SIGPROC_INPUT_HPP

#include <string>
#include "aa_sigproc_input.hpp"
#include "aa_py_filterbank_metadata.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {
      aa_sigproc_input* aa_py_sigproc_input(char const* path);
      void aa_py_sigproc_input_del(aa_sigproc_input const*const obj);
      aa_py_filterbank_metadata_struct aa_py_sigproc_input_read_metadata(aa_sigproc_input *const obj);
      bool aa_py_sigproc_input_read_signal(aa_sigproc_input *const obj);
      unsigned short const* aa_py_sigproc_input_input_buffer_modifiable(aa_sigproc_input const*const obj);
    } // extern "C"
  } // namespace python
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PY_SIGPROC_INPUT_HPP
