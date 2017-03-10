//Added by Karel Adamek

#ifndef __BLN__
#define __BLN__

extern void BLN_init(void);
extern int BLN(float *d_input, float *d_MSD, int CellDim_x, int CellDim_y, int nDMs, int nTimesamples, int offset, float multiplier);

#endif
