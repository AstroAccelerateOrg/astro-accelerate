#ifndef SKA_ASTROACCELERATE_GETDATA_H_
#define SKA_ASTROACCELERATE_GETDATA_H_

void get_recorded_data(FILE **fp, int nsamp, int nchans, int nbits, unsigned short **input_buffer, size_t *inputsize);

#endif
