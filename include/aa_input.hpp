#ifndef ASTRO_ACCELERATE_AA_INPUT_HPP
#define ASTRO_ACCELERATE_AA_INPUT_HPP

#include "aa_filterbank_metadata.hpp"

namespace astroaccelerate {

  /**
   * \class aa_input aa_input.hpp "include/aa_input.hpp"
   * \brief Base class for implementing input data source, which must be implemented in derived classes.
   * \details The base class virtual implmentations return trivial and empty objects.
   * \author Cees Carels.
   * \date 5 November 2018.
   */
  class aa_input {
  public:

    /** \brief Constructor for aa_input.*/    
    aa_input() {

    }
    
    /** \brief Destructor for aa_input. */
    ~aa_input() {

    }
    
    /** Method to open input source for reading. */
    virtual bool open() {
      return false;
    }
    
    /** Method to close input source after reading. */
    virtual bool close() {
      return false;
    }
    
    /** Method to read data from input source. */
    virtual aa_filterbank_metadata read_metadata() {
      if(!isopen) {
	open();
      }
      aa_filterbank_metadata empty;
      return empty;
    }
  protected:
    bool isopen; /** Flag to indicate whether the input source is open. */
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_INPUT_HPP
