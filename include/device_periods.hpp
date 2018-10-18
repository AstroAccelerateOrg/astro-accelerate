#ifndef ASTRO_ACCELERATE_DEVICE_PERIODS_HPP
#define ASTRO_ACCELERATE_DEVICE_PERIODS_HPP

#include "gpu_timer.hpp"

extern void GPU_periodicity(int range, int nsamp, int max_ndms, int processed, float sigma_cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp, int nHarmonics, int candidate_algorithm, int enable_outlier_rejection, float bln_sigma_constant);

#endif








