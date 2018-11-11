#include "aa_input.hpp"

namespace astroaccelerate {

aa_input::aa_input() {

}

aa_input::~aa_input() {

}

bool aa_input::open() {
  return false;
}

bool aa_input::close() {
  return false;
}

aa_filterbank_metadata aa_input::read_metadata() {
    if(!isopen) {
        open();
    }
    aa_filterbank_metadata empty;
    return empty;
}

} //namespace astroaccelerate
