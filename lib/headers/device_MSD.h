#include "headers/device_MSD_Configuration.h"

#ifndef ASTROACCELERATE_MSD_H_
#define ASTROACCELERATE_MSD_H_

extern void MSD_init(void);
extern int MSD_normal(float *d_input, float *d_MSD, float *d_temp, MSD_Configuration *MSD_conf);
extern int MSD_normal(float *d_input, float *d_MSD, int nDMs, int nTimesamples, int offset);
extern int MSD_normal_continuous(float *d_input, float *d_MSD, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf);
extern int MSD_normal_continuous(float *d_input, float *d_MSD, float *d_previous_partials, int nDMs, int nTimesamples, int offset);
extern int MSD_outlier_rejection(float *d_input, float *d_MSD, float *d_temp, MSD_Configuration *MSD_conf, float sigma_outlier_rejection_multiplier);
extern int MSD_outlier_rejection(float *d_input, float *d_MSD, int nDMs, int nTimesamples, int offset, float sigma_outlier_rejection_multiplier);
extern int MSD_outlier_rejection_continuous(float *d_input, float *d_MSD, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float sigma_outlier_rejection_multiplier);
extern int MSD_outlier_rejection_continuous(float *d_input, float *d_MSD, float *d_previous_partials, int nDMs, int nTimesamples, int offset, float sigma_outlier_rejection_multiplier);
extern int MSD_outlier_rejection_grid(float *d_input, float *d_MSD, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float sigma_outlier_rejection_multiplier);
extern int MSD_grid_outlier_rejection(float *d_input, float *d_MSD, int CellDim_x, int CellDim_y, int nDMs, int nTimesamples, int offset, float multiplier);

#endif
