#ifndef ASTRO_ACCELERATE_AA_INPUT_HPP
#define ASTRO_ACCELERATE_AA_INPUT_HPP

#include "aa_filterbank_metadata.hpp"

/**
 * aa_input provides a uniform interface to the data stream, which must be implemented in derived classes.
 */

namespace astroaccelerate {

class aa_input {
public:
    aa_input();
    ~aa_input();
    
    virtual bool open();
    virtual bool close();
    virtual aa_filterbank_metadata read_metadata();
protected:
    bool isopen;
};

} //namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_INPUT_HPP
