#ifndef ASTRO_ACCELERATE_AA_SIGPROC_INPUT_HPP
#define ASTRO_ACCELERATE_AA_SIGPROC_INPUT_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "aa_input.hpp"

namespace astroaccelerate {

  /**
   * \class aa_sigproc_input aa_sigproc_input.hpp "include/aa_sigproc_input.hpp"
   * \brief The aa_sigproc_input class is used to parse a sigproc (.fil) file.
   * \details The aa_sigproc_input class is used to parse a sigproc (.fil) file.
   * \details Implementation details provided in source src/aa_sigproc_input.cpp.
   * \details Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs"
   * \details Source: http://sigproc.sourceforge.net/sigproc.pdf
   * \author Cees Carels.
   * \date 5 November 2018.
   */
  class aa_sigproc_input : public aa_input {
  public:

    /** \brief Constructor for aa_sigproc_input, requires a file path. */
    aa_sigproc_input(const std::string &path);

    /** \brief Destructor for aa_sigproc_input. */
    virtual ~aa_sigproc_input();

    /** \brief Method to open the input file. */
    bool open();

    /** \brief Method to close the input file. */
    bool close();

    /**
     * \brief Method to read only the metadata from the input file.
     * \returns an aa_filterbank_metadata object with the metadata from the input file. */
    aa_filterbank_metadata read_metadata();

    /** \brief Method to read only the input data from the input file. */
    bool read_signal();

    /** \returns The input data from the telescope. */
    const std::vector<unsigned short>& input_buffer() const {
      return m_input_buffer;
    }

    /** \returns A modifiable reference to the input data from the telescope. */
    std::vector<unsigned short>& input_buffer_modifiable() {
      return m_input_buffer;
    }

    /** \brief Method to check if the input data from the input file has already been read.
     * \returns A boolean flag to indicate if the data has been read (true) or not (false). */
    bool did_read_signal() const {
      return data_is_read;
    }
    
    
  private:
    bool get_file_data(aa_filterbank_metadata &metadata); /** Reads the metadata from the filterbank input file. */
    
    std::vector<unsigned short> m_input_buffer; /** Stores the data in the sigproc file. */
    
    template <typename T>
    bool get_recorded_data(std::vector<T> &input_buffer);
    bool header_is_read; /** Flag to indicate whether the input file header (metadata) has been read. */
    bool data_is_read; /** Flag to indicate whether the input file data has been read. */
    
    aa_filterbank_metadata m_meta; /** Stores the metadata associated with the input file. */
    
    std::string file_path; /** Stores the file path to the input file. */
    FILE *fp; /** File pointer to the input file. */
    
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_SIGPROC_INPUT_HPP
