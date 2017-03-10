#ifndef ASTROACCELERATE_THRESHOLD_H_
#define ASTROACCELERATE_THRESHOLD_H_

#include<vector>
#include "device_BC_plan.h"


extern void THR_init(void);
extern int THRESHOLD(float *d_input, ushort *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nDMs, int nTimesamples, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int max_list_size);

#endif

