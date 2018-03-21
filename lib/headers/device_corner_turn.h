#ifndef ASTROACCELERATE_CORNERTURN_H_
#define ASTROACCELERATE_CORNERTURN_H_

extern void corner_turn(unsigned short *d_input, float *d_output, int nchans, int nsamp, cudaStream_t stream);
extern void corner_turn(float *d_input, float *d_output, int primary_size, int secondary_size);
extern void corner_turn_SM(float *d_input, float *d_output, int primary_size, int secondary_size);


#endif

