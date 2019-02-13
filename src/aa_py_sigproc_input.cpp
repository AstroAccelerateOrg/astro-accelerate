#include "aa_py_sigproc_input.hpp"
#include <iostream>
namespace astroaccelerate {
  namespace python {
    extern "C" {
      aa_sigproc_input* aa_py_sigproc_input(char const* path) {
	return new aa_sigproc_input(path);
      }

      void aa_py_sigproc_input_del(aa_sigproc_input const*const obj) {
	delete obj;
      }

      aa_py_filterbank_metadata_struct aa_py_sigproc_input_read_metadata(aa_sigproc_input *const obj) {
	aa_filterbank_metadata meta = obj->read_metadata();
	aa_py_filterbank_metadata_struct meta_obj;
	meta_obj.m_tstart = meta.tstart();
	meta_obj.m_tsamp = meta.tsamp();
	meta_obj.m_fch1 = meta.fch1();
	meta_obj.m_foff = meta.foff();
	
	meta_obj.m_nbits = meta.nbits();
	meta_obj.m_nsamples = meta.nsamples();
	meta_obj.m_nchans = meta.nchans();
	return meta_obj;
      }

      bool aa_py_sigproc_input_read_signal(aa_sigproc_input *const obj) {
	return obj->read_signal();
      }

      unsigned short const* aa_py_sigproc_input_input_buffer(aa_sigproc_input const*const obj) {
	return obj->input_buffer().data();
      }
    }
  }
}
