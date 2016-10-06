#ifndef __MSD_LIMITED__
#define __MSD_LIMITED__

extern void MSD_limited_init(void);
extern int MSD_limited(float *d_input, float *d_MSD, int nDMs, int nTimesamples, int offset);

#endif
