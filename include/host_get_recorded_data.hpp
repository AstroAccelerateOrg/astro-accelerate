#ifndef ASTRO_ACCELERATE_HOST_GET_RECRODED_DATA_HPP
#define ASTRO_ACCELERATE_HOST_GET_RECORDED_DATA_HPP

void get_recorded_data(FILE**           fp,
                       int              nsamp,
                       int              nchans,
                       int              nbits,
                       unsigned short** input_buffer,
                       size_t*          inputsize);

#endif
