#ifndef ASTRO_ACCELERATE_AA_DEVICE_SPS_LONG_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SPS_LONG_HPP

#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "aa_params.hpp"
#include "aa_device_BC_plan.hpp"
#include "aa_device_SPS_long_kernel.hpp"

namespace astroaccelerate {
  
  /** \brief Function that calculates the amount of memory required to perform the SPS (analysis). */
  size_t Get_memory_requirement_of_SPS();
  extern void PD_SEARCH_LONG_init();
  extern int SPDT_search_long_MSD_plane(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD_interpolated, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nTimesamples, int nDMs);	

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SPS_LONG_HPP
