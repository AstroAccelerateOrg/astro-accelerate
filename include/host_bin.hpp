#ifndef ASTRO_ACCELERATE_HOST_BIN_HPP
#define ASTRO_ACCELERATE_HOST_BIN_HPP

void bin(size_t binsize,
         float *bin_buffer,
         size_t inputsize,
         float *input_buffer,
         int    nchans,
         int    nsamps);

#endif
