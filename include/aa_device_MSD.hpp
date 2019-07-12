#ifndef ASTRO_ACCELERATE_AA_DEVICE_MSD_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_MSD_HPP

#include <vector>

#include "aa_device_MSD_Configuration.hpp"
#include "aa_device_MSD_shared_kernel_functions.hpp"
#include "aa_device_MSD_normal_kernel.hpp"
#include "aa_device_MSD_outlier_rejection_kernel.hpp"

namespace astroaccelerate {

extern void MSD_init(void);
extern int MSD_normal(float *d_MSD, float *d_input, float *d_temp, MSD_Configuration *MSD_conf);
extern int MSD_normal(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset);
extern int MSD_normal_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf);
extern int MSD_normal_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset);
extern int MSD_outlier_rejection(float *d_MSD, float *d_input, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier);
extern int MSD_outlier_rejection(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier);
extern int MSD_outlier_rejection_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier);
extern int MSD_outlier_rejection_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier);
extern int MSD_outlier_rejection_grid(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier);
extern int MSD_grid_outlier_rejection(float *d_MSD, float *d_input, int CellDim_x, int CellDim_y, int nTimesamples, int nDMs, int offset, float multiplier);


extern void Find_MSD(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection);
extern void Find_MSD(float *d_MSD, float *d_input, float *d_MSD_workarea, MSD_Configuration *conf, float OR_sigma_multiplier, int enable_outlier_rejection);
extern void Find_MSD_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection);
extern void Find_MSD_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_MSD_workarea, MSD_Configuration *conf, float OR_sigma_multiplier, int enable_outlier_rejection);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_MSD_HPP
