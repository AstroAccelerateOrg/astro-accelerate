#ifndef ASTRO_ACCELERATE_AA_DEVICE_THRESHOLD_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_THRESHOLD_HPP

#include <vector>
#include "aa_device_BC_plan.hpp"

namespace astroaccelerate {

  extern void THR_init(void);
  extern int SPDT_threshold(float *d_input, ushort *d_input_taps, unsigned int *d_output_list_DM, unsigned int *d_output_list_TS, float *d_output_list_SNR, unsigned int *d_output_list_BW, int *gmem_pos, float threshold, int nDMs, int nTimesamples, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int max_list_size);

  extern int Threshold_for_periodicity_old(float *d_input, ushort *d_input_harms, float *d_output_list, int *gmem_pos, float *d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int inBin, int max_list_size);

  extern int Threshold_for_periodicity(float *d_input, ushort *d_input_harms, float *d_output_list, int *gmem_pos, float *d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int inBin, int max_list_size);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_THRESHOLD_HPP

