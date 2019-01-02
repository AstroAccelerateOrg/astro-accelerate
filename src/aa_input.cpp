#include "aa_input.hpp"

namespace astroaccelerate {

  /**
   * Constructor for aa_input.
   */
  aa_input::aa_input() {

  }
  
  /**
   * Destructor for aa_input.
   */
  aa_input::~aa_input() {

  }

  /**
   * Method to open the input for reading.
   */
  bool aa_input::open() {
    return false;
  }

  /**
   * Method to close the input after reading.
   */
  bool aa_input::close() {
    return false;
  }

  /**
   * Method to read the data from the input.
   */
  aa_filterbank_metadata aa_input::read_metadata() {
    if(!isopen) {
      open();
    }
    aa_filterbank_metadata empty;
    return empty;
  }

} //namespace astroaccelerate
