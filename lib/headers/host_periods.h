#ifndef ASTROACCELERATE_PERIODS_H_
#define ASTROACCELERATE_PERIODS_H_

void periodicity(int range, int nsamp, int max_ndms, int processed, int num_boots, int num_trial_bins, int navdms, float narrow, float wide, int nsearch, float aggression, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp);

#endif

