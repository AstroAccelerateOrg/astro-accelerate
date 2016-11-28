#ifndef __THRESHOLD__
#define __THRESHOLD__

extern void THR_init(void);
extern int THRESHOLD(float *d_input, unsigned char *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nDMs, int nTimesamples, int offset, int max_list_size);
extern int THRESHOLD_ignore(float *d_input, unsigned char *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nTaps, int nDMs, int nTimesamples, int offset, int max_list_size);

#endif

