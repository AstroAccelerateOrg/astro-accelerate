#ifndef SKA_ASTROACCELERATE_MSD_LIMITED_H_
#define SKA_ASTROACCELERATE_MSD_LIMITED_H_

extern void MSD_limited_init(void);
extern int MSD_limited(float *d_input, float *d_MSD, int nDMs, int nTimesamples, int offset);

#endif
