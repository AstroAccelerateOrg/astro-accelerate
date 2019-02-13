#include "aa_py_filterbank_metadata.hpp"
#include "aa_filterbank_metadata.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {
      aa_py_filterbank_metadata_struct* aa_py_filterbank_metadata(const double tstart,
								  const double tsamp,
								  const int nbits,
								  const int nsamples,
								  const double fch1,
								  const double foff,
								  const int nchans) {
	aa_py_filterbank_metadata_struct* obj = new aa_py_filterbank_metadata_struct;

	obj->m_tstart = tstart;
	obj->m_tsamp = tsamp;
	obj->m_fch1 = fch1;
	obj->m_foff = foff;

	obj->m_nbits = nbits;
	obj->m_nsamples = nsamples;
	obj->m_nchans = nchans;

	return obj;
      }
    
      void aa_py_filterbank_metadata_delete(aa_py_filterbank_metadata_struct const*const obj) {
	delete obj;
      }
    
      double aa_py_filterbank_metadata_tstart(aa_py_filterbank_metadata_struct const*const obj) {
	return obj->m_tstart;
      }
    } // extern "C"
  } // namespace python
} // namespace astroaccelerate
