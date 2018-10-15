#ifndef ASTRO_ACCELERATE_DEVICE_LOAD_DATA_HPP
#define ASTRO_ACCELERATE_DEVICE_LOAD_DATA_HPP

void load_data(int             i,
               int*            inBin,
               unsigned short* device_pointer,
               unsigned short* host_pointer,
               int             t_processed,
               int             maxshift,
               int             nchans,
               float*          dmshifts);

#endif
