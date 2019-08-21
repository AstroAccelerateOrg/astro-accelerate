// James Sharpe's peak finding code

#ifndef ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_HPP

#include <vector>
#include <npp.h>

#include "aa_params.hpp"
#include "aa_device_peak_find_kernel.hpp"
#include "aa_device_BC_plan.hpp"

namespace astroaccelerate {

  extern void SPDT_peak_find(float *d_output_SNR, ushort *d_output_taps, unsigned int *d_peak_list_DM, unsigned int *d_peak_list_TS, float *d_peak_list_SNR, unsigned int *d_peak_list_BW, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration);

  extern void PEAK_FIND_FOR_FDAS(float *d_ffdot_plane, float *d_peak_list, float *d_MSD, int nDMs, int nTimesamples, float threshold, unsigned int max_peak_size, unsigned int *gmem_peak_pos, float DM_trial);

  extern void Peak_find_for_periodicity_search_old(float *d_input_SNR, ushort *d_input_harmonics, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, float *d_MSD, int DM_shift, int inBin);

  extern void Peak_find_for_periodicity_search(float *d_input_SNR, ushort *d_input_harmonics, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, float* d_MSD, int DM_shift, int inBin);

  extern void SPDT_peak_find_stencil_7x7(float *d_output_SNR, ushort *d_output_taps, unsigned int *d_peak_list_DM, unsigned int *d_peak_list_TS, float *d_peak_list_SNR, unsigned int *d_peak_list_BW, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration);

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_HPP
