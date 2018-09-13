#ifndef ASTROACCELERATE_MSD_GRID_H_
#define ASTROACCELERATE_MSD_GRID_H_

extern void MSD_grid_init(void);
extern int MSD_grid(float *d_input, float *d_output, int CellDim_x, int CellDim_y, int nDMs, int nTimesamples);

#endif

