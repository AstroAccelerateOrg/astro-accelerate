#include <vector>
#include "headers/device_MSD_Configuration.h"

#ifndef ASTROACCELERATE_MSD_H_
#define ASTROACCELERATE_MSD_H_

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

extern void Get_MSD_plane_profile_memory_requirements(size_t *MSD_profile_size_in_bytes, size_t *MSD_DIT_profile_size_in_bytes, size_t *workarea_size_in_bytes, size_t primary_dimension, size_t secondary_dimension, std::vector<int> *boxcar_widths);

extern void MSD_plane_profile(float *d_MSD_interpolated, float *d_input_data, float *d_MSD_DIT_previous, float *workarea, bool high_memory, size_t primary_dimension, size_t secondary_dimension, std::vector<int> *boxcar_widths, float tstart, float dm_low, float dm_high, float OR_sigma_multiplier, int enable_outlier_rejection, bool perform_continuous);

extern void MSD_plane_profile_boxcars(float *d_input_data, size_t nTimesamples, size_t nDMs, std::vector<int> *boxcar_widths, float OR_sigma_multiplier, float dm_low, float dm_high, float tstart);
#endif
