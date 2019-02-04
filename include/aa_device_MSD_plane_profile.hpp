#ifndef ASTRO_ACCELERATE_AA_DEVICE_MSD_PLANE_PROFILE_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_MSD_PLANE_PROFILE_HPP

#include <vector>

namespace astroaccelerate {

extern void MSD_plane_profile(float *d_MSD_interpolated, float *d_input_data, float *d_MSD_DIT_previous, float *workarea, bool high_memory, size_t primary_dimension, size_t secondary_dimension, std::vector<int> *boxcar_widths, float tstart, float dm_low, float dm_high, float OR_sigma_multiplier, int enable_outlier_rejection, bool perform_continuous, double *total_time, double *dit_time, double *MSD_time);
extern void Get_MSD_plane_profile_memory_requirements(size_t *MSD_profile_size_in_bytes, size_t *MSD_DIT_profile_size_in_bytes, size_t *workarea_size_in_bytes, size_t primary_dimension, size_t secondary_dimension, std::vector<int> *boxcar_widths);
extern void MSD_plane_profile_boxcars(float *d_input_data, size_t nTimesamples, size_t nDMs, std::vector<int> *boxcar_widths, float OR_sigma_multiplier, float dm_low, float dm_high, float tstart);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_MSD_PLANE_PROFILE_HPP

