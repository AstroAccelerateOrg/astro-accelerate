#ifndef ASTRO_ACCELERATE_HOST_MSD_STRATEGY_HPP
#define ASTRO_ACCELERATE_HOST_MSD_STRATEGY_HPP

typedef struct {
  unsigned long int maxtimesamples;
} MSD_info;

void stratagy_MSD(int                ndms,
                  float              max_boxcar_width_in_sec,
                  float              tsamp,
                  int                nTimesamples,
                  unsigned long int *info,
                  size_t *           MSD_profile_size_in_bytes,
                  int *              MSD_nDIT_widths);

#endif
