#ifndef ASTRO_ACCELERATE_DEVICE_MSD_GRID_HPP
#define ASTRO_ACCELERATE_DEVICE_MSD_GRID_HPP

//It looks as though these function prototypes are declared but not implemented or used anywhere
extern void MSD_grid_init(void);
extern int MSD_grid(float *d_input, float *d_output, int CellDim_x, int CellDim_y, int nDMs, int nTimesamples);

#endif

