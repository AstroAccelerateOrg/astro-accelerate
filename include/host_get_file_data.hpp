#ifndef ASTRO_ACCELERATE_HOST_GET_FILE_DATA_HPP
#define ASTRO_ACCELERATE_HOST_GET_FILE_DATA_HPP

void get_file_data(FILE** fp,
                   int*   nchans,
                   int*   nsamples,
                   int*   nsamp,
                   int*   nifs,
                   int*   nbits,
                   float* tsamp,
                   float* tstart,
                   float* fch1,
                   float* foff);

#endif
