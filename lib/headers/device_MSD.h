#include <vector>
#include "headers/device_MSD_Configuration.h"

#ifndef ASTROACCELERATE_MSD_H_
#define ASTROACCELERATE_MSD_H_

extern void MSD_init(void);
extern int MSD_normal(float *d_MSD, float *d_input, float *d_temp, MSD_Configuration *MSD_conf, cudaStream_t streams);
extern int MSD_normal(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, cudaStream_t streams);
extern int MSD_normal_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, cudaStream_t streams);
extern int MSD_normal_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, cudaStream_t streams);
extern int MSD_outlier_rejection(float *d_MSD, float *d_input, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier, cudaStream_t streams);
extern int MSD_outlier_rejection(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, cudaStream_t streams);
extern int MSD_outlier_rejection_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier, cudaStream_t streams);
extern int MSD_outlier_rejection_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, cudaStream_t streams);
extern int MSD_outlier_rejection_grid(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier, cudaStream_t streams);
extern int MSD_grid_outlier_rejection(float *d_MSD, float *d_input, int CellDim_x, int CellDim_y, int nTimesamples, int nDMs, int offset, float multiplier, cudaStream_t streams);


extern void Find_MSD(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection, cudaStream_t streams);
extern void Find_MSD(float *d_MSD, float *d_input, float *d_MSD_workarea, MSD_Configuration *conf, float OR_sigma_multiplier, int enable_outlier_rejection, cudaStream_t streams);
extern void Find_MSD_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection, cudaStream_t streams);
extern void Find_MSD_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_MSD_workarea, MSD_Configuration *conf, float OR_sigma_multiplier, int enable_outlier_rejection, cudaStream_t streams);
#endif
