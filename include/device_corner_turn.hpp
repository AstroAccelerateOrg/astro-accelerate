#ifndef ASTRO_ACCELERATE_DEVICE_CORNER_TURN_HPP
#define ASTRO_ACCELERATE_DEVICE_CORNER_TURN_HPP

extern void
            corner_turn(unsigned short *d_input, float *d_output, int nchans, int nsamp);
extern void corner_turn(float *d_input,
                        float *d_output,
                        int    primary_size,
                        int    secondary_size);
extern void corner_turn_SM(float *d_input,
                           float *d_output,
                           int    primary_size,
                           int    secondary_size);

#endif
