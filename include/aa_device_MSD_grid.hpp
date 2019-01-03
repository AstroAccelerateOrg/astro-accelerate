#ifndef ASTRO_ACCELERATE_AA_DEVICE_MSD_GRID_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_MSD_GRID_HPP

namespace astroaccelerate {

//It looks as though these function prototypes are declared but not implemented or used anywhere
extern void MSD_grid_init(void);
extern int MSD_grid(float *d_input, float *d_output, int CellDim_x, int CellDim_y, int nDMs, int nTimesamples);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_MSD_GRID_HPP
