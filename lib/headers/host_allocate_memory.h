#ifndef ASTROACCELERATE_ALLOCATEMEMORY_H_
#define ASTROACCELERATE_ALLOCATEMEMORY_H_

#include "headers/device_DDTR_Plan.h"

void allocate_memory_cpu_input(unsigned short **input_buffer, size_t nsamp, size_t nchans);

void allocate_memory_cpu_output(float ****output_buffer, DDTR_Plan *DDTR_plan);

void allocate_memory_gpu(unsigned short **d_input, float **d_output, DDTR_Plan *DDTR_plan);

#endif
