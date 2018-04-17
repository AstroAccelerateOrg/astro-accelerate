#ifndef ASTROACCELERATE_LOAD_DATA_H_
#define ASTROACCELERATE_LOAD_DATA_H_

void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short *host_pointer, unsigned short *host_pointer_pinned, int t_processed, int maxshift, int nchans, float *dmshifts, cudaStream_t stream);

#endif

