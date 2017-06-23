//Added by Karel Adamek

#ifndef __MSD_BLN_GRID__
#define __MSD_BLN_GRID__

extern void MSD_BLN_grid_init(void);
extern int MSD_BLN_grid(float *d_input, float *d_MSD, int CellDim_x, int CellDim_y, int nDMs, int nTimesamples, int offset, float multiplier);

#endif
