#ifndef ASTRO_ACCELERATE_AA_SIGPROC_INPUT_HPP
#define ASTRO_ACCELERATE_AA_SIGPROC_INPUT_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "aa_input.hpp"
#include "aa_host_malloc.hpp"

/**
 * The aa_sigproc_input class is used to parse a sigproc (.fil) file.
 * Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs"
 * Source: http://sigproc.sourceforge.net/sigproc.pdf
 */

class aa_sigproc_input : public aa_input {
public:
    aa_sigproc_input(const std::string &path);
    ~aa_sigproc_input();
    bool open();
    bool close();
    aa_filterbank_metadata read_metadata();
    bool read_telescope();
    const std::vector<unsigned short>& input_buffer() const {
        return m_input_buffer;
    }
    
    
private:
    bool get_file_data(aa_filterbank_metadata &metadata);
    
    std::vector<unsigned short> m_input_buffer;
    
    template <typename T>
    bool get_recorded_data(std::vector<T> &input_buffer);
    bool header_is_read;
    
    aa_filterbank_metadata m_meta;
    aa_malloc              m_mem;
    
    std::string file_path;
    FILE *fp;
    
};

#endif // ASTRO_ACCELERATE_AA_SIGPROC_INPUT_HPP
