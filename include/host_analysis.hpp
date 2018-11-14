#ifndef ASTRO_ACCELERATE_HOST_ANALYSIS_HPP
#define ASTRO_ACCELERATE_HOST_ANALYSIS_HPP

namespace astroaccelerate {

void analysis_CPU(int i, float tstart, int t_processed, int max_ndms, int *ndms, int *outBin, float cutoff, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, float max_boxcar_width_in_sec);

} //namespace astroaccelerate

#endif

