#include "aa_sigproc_input.hpp"

#include "aa_log.hpp"

namespace astroaccelerate {

  /** \brief Constructor for aa_sigproc_input. */
  aa_sigproc_input::aa_sigproc_input(const std::string &path) : header_is_read(false), data_is_read(false), file_path(path) {
    isopen = false;
  }

  /** \brief Destructor for aa_sigproc_input. */
  aa_sigproc_input::~aa_sigproc_input() {
    if(isopen) {
      close();
    }
  }

  /**
   * \brief Method to open the sigproc input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::open() {
    fp = fopen(file_path.c_str(), "rb");
    if(fp == NULL) {
      return false;
    }
    isopen = true;
    return true;
  }

  /**
   * \brief Closes the input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::close() {
    if(isopen) {
      if(fclose(fp) == -1) {
	return false;
      }
      else {
	isopen = false;
      }
    }
    return true;
  }

  /**
   * \brief Method to read the metadata from the sigproc input file.
   * \returns An aa_filterbank_metadata object containing the metadata read from the sigproc input file. If the data could not be read, a trivial instance is returned.
   */
  aa_filterbank_metadata aa_sigproc_input::read_metadata() {
    if(!isopen) {
      if(!open()) {
	aa_filterbank_metadata empty;
	return empty;
      }
    }
    
    aa_filterbank_metadata metadata;
    get_file_data(metadata);
    header_is_read = true;
    return metadata;
  }

  /**
   * \brief If the file is open, and the header has been read, and the data has not yet been read, then read the input data from the input file.
   * \details Reading the telescope input data can only be performed once, after which this method will always return false.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   * \warning The method will return true only once, that is the first time the data are read from the input successfully. At this point the input_buffer should be checked for data. 
   */
  bool aa_sigproc_input::read_signal() {
    if(!isopen || !header_is_read || data_is_read) {
      return false;
    }
    
    get_recorded_data(m_input_buffer);
    return true;
  }

  /**
   * \brief Method to read the input data from the sigproc input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::get_file_data(aa_filterbank_metadata &metadata) {
    double az_start = 0;
    double za_start = 0;
    double src_raj = 0;
    double src_dej = 0;
    double tstart = 0;
    double tsamp = 0;
    double refdm = 0;
    double period = 0;
    double fch1 = 0;
    double foff = 0;
    double fchannel = 0;
    
    int telescope_id = 0;
    int machine_id = 0;
    int data_type = 0;
    int barycentric = 0;
    int pulsarcentric = 0;
    int nbits = 0;
    int nsamples = 0;
    int nchans = 0;
    int nifs = 0;
    
    char *string = (char *) malloc(80 * sizeof(char));
    
    int nchar = 0;
    int nbytes = sizeof(int);
    
    while (1) {
      strcpy(string, "ERROR");
      if (fread(&nchar, sizeof(int), 1, fp) != 1) {
	fprintf(stderr, "\nError while reading file\n");
	return false;
      }
      if (feof(fp)) {
	return false;
      }
        
      if (nchar > 1 && nchar < 80) {
	if (fread(string, nchar, 1, fp) != 1) {
	  fprintf(stderr, "\nError while reading file\n");
	  return false;
	}
            
	string[nchar] = '\0';
	// For debugging only
	//printf("From .fil header:\t%d\t%s\n", nchar, string);
	nbytes += nchar;
            
	if (strcmp(string, "HEADER_END") == 0) {
	  break;
	}
            
	if(strcmp(string, "telescope_id") == 0) {
	  if (fread(&telescope_id, sizeof(telescope_id), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "machine_id") == 0) {
	  if (fread(&machine_id, sizeof(machine_id), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "data_type") == 0) {
	  if (fread(&data_type, sizeof(data_type), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "barycentric") == 0) {
	  printf("barycentric\n");
	  if (fread(&barycentric, sizeof(barycentric), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "pulsarcentric") == 0) {
	  if (fread(&pulsarcentric, sizeof(pulsarcentric), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "az_start") == 0) {
	  if (fread(&az_start, sizeof(az_start), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "za_start") == 0) {
	  if (fread(&za_start, sizeof(za_start), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "src_raj") == 0) {
	  if (fread(&src_raj, sizeof(src_raj), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "src_dej") == 0) {
	  if (fread(&src_dej, sizeof(src_dej), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "tstart") == 0) {
	  if (fread(&tstart, sizeof(tstart), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "tsamp") == 0) {
	  if (fread(&tsamp, sizeof(tsamp), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "nbits") == 0) {
	  if (fread(&nbits, sizeof(nbits), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "nsamples") == 0) {
	  if (fread(&nsamples, sizeof(nsamples), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "fch1") == 0) {
	  if (fread(&fch1, sizeof(fch1), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "foff") == 0) {
	  if (fread(&foff, sizeof(foff), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "fchannel") == 0) {
	  if (fread(&fchannel, sizeof(fchannel), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "nchans") == 0) {
	  if (fread(&nchans, sizeof(nchans), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "nifs") == 0) {
	  if (fread(&nifs, sizeof(nifs), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "refdm") == 0) {
	  if (fread(&refdm, sizeof(refdm), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "period") == 0) {
	  if (fread(&period, sizeof(period), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
      }
    }
    
    free(string);
    
    // Check that we are working with one IF channel
    if (nifs != 1) {
      LOG(log_level::error, "Astro-Accelerate can only work with one IF channel.");
      return false;
    }
    
    fpos_t file_loc;
    //Position of end of header
    if (fgetpos(fp, &file_loc) != 0){
	    LOG(log_level::error, "Could not get position of the header end.");
	    return false;
    }
    
    unsigned long int nsamp = 0;

    // Getting number of time-samples based on file size
	long int data_start = ftell(fp);
	if (data_start == -1){
		LOG(log_level::error, "Wrong return of the start position indicator of the file (stream).");
		return false;
	}


    if (fseek(fp, 0, SEEK_END) != 0) {
      LOG(log_level::error, "Failed to seek to the end of the data file");
      return false;
    }

    long int exp_total_data = ftell(fp);
	if (exp_total_data == -1){
		LOG(log_level::error, "Wrong return of the end position indicator of the file (stream).");
		return false;
	}

    exp_total_data = exp_total_data - data_start;

    if (fseek(fp, data_start, SEEK_SET) != 0){
	    LOG(log_level::error, "Failed to seek to the data_start.");
	    return false;
    }

    if (( nbits ) == 32) {
      nsamp = exp_total_data/((nchans)*4);
    }
    else if (( nbits ) == 16) {
      nsamp = exp_total_data/((nchans)*2);
    }
    else if (( nbits ) == 8) {
      nsamp = exp_total_data/((nchans));
    }
    else if (( nbits ) == 4) {
      nsamp = 2*exp_total_data/((nchans));
    }
    else if (( nbits ) == 2) {
      nsamp = 4*exp_total_data/((nchans));
    }
    else if (( nbits ) == 1) {
      nsamp = 8*exp_total_data/((nchans));
    }
    else {
      LOG(log_level::error, "Currently this code only runs with 1, 2, 4, 8 and 16 bit data.");
      return false;
    }
    
    // Move the file pointer back to the end of the header
    if (fsetpos(fp, &file_loc) != 0){
	LOG(log_level::error, "Could not set position to the end of the header in the file.");
	return false;
    }
    
    aa_filterbank_metadata meta(telescope_id,
                                machine_id,
                                data_type,
                                "unknown",
                                "unknown",
                                barycentric,
                                pulsarcentric,
                                az_start,
                                za_start,
                                src_raj,
                                src_dej,
                                tstart,
                                tsamp,
                                nbits,
                                nsamp,
                                fch1,
                                foff,
                                0,
                                fchannel,
                                0,
                                nchans,
                                nifs,
                                refdm,
                                period);
    metadata = meta;
    m_meta   = meta;
    
    return true;
  }

  template <typename T>
  bool aa_sigproc_input::get_recorded_data(std::vector<T> &input_buffer) {
    const size_t inputsize = (size_t)m_meta.nsamples() * (size_t)m_meta.nchans();
    input_buffer.resize(inputsize);
    int c;
    
    unsigned long int total_data;
    const int nchans = m_meta.nchans();
    const int nbits = m_meta.nbits();
    //{{{ Load in the raw data from the input file and transpose
    if (nbits == 32) {
      // Allocate a tempory buffer to store a line of frequency data
      float *temp_buffer = (float *) malloc(nchans * sizeof(float));
        
      // Allocate floats to hold the mimimum and maximum value in the input data
      float max = -10000000000;
      float min = 10000000000;
        
      // Allocate a variable to hold the file pointer position
      fpos_t file_loc;
        
      // Set the file pointer position
      fgetpos(fp, &file_loc);
        
      // Find the minimum and maximum values in the input file.
      while (!feof(fp)) {
	if (fread(temp_buffer, sizeof(float), nchans, fp) != (size_t)nchans)
	  break;
	for (c = 0; c < nchans; c++) {
	  if(temp_buffer[c] > max) max = temp_buffer[c];
	  if(temp_buffer[c] < min) min = temp_buffer[c];
	}
      }
        
      // Calculate the bin size in a distribution of unsigned shorts for the input data.
      float bin = (max - min) / 65535.0f;
        
      printf("\n Conversion Parameters: %f\t%f\t%f", min, max, bin);
        
      // Move the file pointer back to the end of the header
      fsetpos(fp, &file_loc);
        
      // Read in the data, scale to an unsigned short range, transpose it and store it in the input buffer
      total_data = 0;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(float), nchans, fp) != (size_t)nchans)
	  break;
	for (c = 0; c < nchans; c++) {
	  ( input_buffer )[c + total_data * ( nchans )] = (unsigned short) ((temp_buffer[c]-min)/bin);
	}
	total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else if (nbits == 16) {
      // Allocate a tempory buffer to store a line of frequency data
      unsigned short *temp_buffer = (unsigned short *) malloc(nchans * sizeof(unsigned short));
        
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned short), nchans, fp) != (size_t)nchans)
	  break;
	for (c = 0; c < nchans; c++) {
	  ( input_buffer )[c + total_data * ( nchans )] = (unsigned short) temp_buffer[c];
	}
	total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else if (nbits == 8) {
        
      // Allocate a tempory buffer to store a line of frequency data
      unsigned char *temp_buffer = (unsigned char *) malloc(nchans * sizeof(unsigned char));
        
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nchans, fp) != (size_t)nchans) {
	  break;
	}
	for (c = 0; c < nchans; c++) {
	  ( input_buffer )[c + total_data * ( nchans )] = (unsigned short) temp_buffer[c];
	}
	total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else if (nbits == 4) {
      // Allocate a temporary buffer to store a line of frequency data
      // each byte stores 2 frequency data
      int nb_bytes = nchans/2;
      unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      // 00001111
      char mask = 0x0f;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
	  break;
	for (c = 0; c < nb_bytes; c++) {
	  // (n >> a) & ( (1 << a) - 1) -> right shift by 'a' bits, then keep the last 'b' bits
	  // Here, we right shift 4 bits and keep the last 4 bits
	  ( input_buffer )[ (c*2) + total_data * ( nchans )]     = (unsigned short)( (temp_buffer[c] >> 4) & mask );
	  // n & ( (1 << a ) - 1)
	  // Here we keep the last 4 bits
	  ( input_buffer )[ (c*2) + 1 + total_data * ( nchans )] = (unsigned short)( temp_buffer[c] & mask );
	}
	total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else if (nbits == 2) {
      // Allocate a temporary buffer to store a line of frequency data
      // each byte stores 4 frequency data
      int nb_bytes = nchans/4;
      unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      // 00001111
      char mask = 0x03;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
	  break;
	for (c = 0; c < nb_bytes; c++) {
	  ( input_buffer )[ (c*4) + total_data * ( nchans )]     = (unsigned short)( (temp_buffer[c] >> 6) & mask );
	  ( input_buffer )[ (c*4) + 1 + total_data * ( nchans )] = (unsigned short)( (temp_buffer[c] >> 4) & mask );
	  ( input_buffer )[ (c*4) + 2 + total_data * ( nchans )] = (unsigned short)( (temp_buffer[c] >> 2) & mask );
	  ( input_buffer )[ (c*4) + 3 + total_data * ( nchans )] = (unsigned short)( temp_buffer[c] & mask );
	}
	total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else if (nbits == 1) {
      // each byte stores 8 frequency data
      int nb_bytes = nchans/8;
        
      // Allocate a temporary buffer to store a line of frequency data
      unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
        
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
        
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
	  break;
            
	for (c = 0; c < nb_bytes; c++) {
	  for(int i=0; i<8; i++) {
	    unsigned char mask =  1 << i;
	    unsigned char masked_char = temp_buffer[c] & mask;
	    unsigned char bit = masked_char >> i;
	    ( input_buffer )[ (c*8) + i + total_data * ( nchans )] = (unsigned short)bit;
	  }
	}
	total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else {
      printf("ERROR: Invalid number of bits in input data.\n");
      return false;
    }

    data_is_read = true;
    return true;
  }

} //namespace astroaccelerate
