// James Sharpe's peak finding code

#ifndef __PEAK_FIND__
#define __PEAK_FIND__

#include<vector>
#include "AstroAccelerate/device_BC_plan.h"

extern void PEAK_FIND(float *d_output_SNR, ushort *d_output_taps, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration);

#endif
