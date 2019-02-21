#ifndef ASTRO_ACCELERATE_AA_PY_FILTERBANK_METADATA_HPP
#define ASTRO_ACCELERATE_AA_PY_FILTERBANK_METADATA_HPP

#include "aa_filterbank_metadata.hpp"

namespace astroaccelerate {
  namespace python {
    extern "C" {

      /**
       * \struct aa_py_filterbank_metadata_struct
       * \brief Wrapper object for aa_filterbank_metadata.
       * \details Please see aa_filterbank_metadata.hpp for library implementation details.
       */
      struct aa_py_filterbank_metadata_struct {
	double m_tstart;
	double m_tsamp;
	double m_fch1;
	double m_foff;

	int m_nbits;
	int m_nsamples;
	int m_nchans;
      };
      
      aa_py_filterbank_metadata_struct* aa_py_filterbank_metadata(const double tstart,
								  const double tsamp,
								  const int nbits,
								  const int nsamples,
								  const double fch1,
								  const double foff,
								  const int nchans);
      
      void aa_py_filterbank_metadata_delete(aa_py_filterbank_metadata_struct const*const obj);
    
      double aa_py_filterbank_metadata_tstart(aa_py_filterbank_metadata_struct const*const obj);
    } // extern "C"
  } // namespace python
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PY_FILTERBANK_METADATA_HPP
