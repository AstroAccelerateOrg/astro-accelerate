#ifndef ASTROACCELERATE_CPUBOOTSTRAP_H_
#define ASTROACCELERATE_CPUBOOTSTRAP_H_

void cpu_blocked_bootstrap(float *data_array, unsigned int num_els, unsigned int num_bins, unsigned int num_boots, float *mean_data, float *cpu_mean, float *cpu_sd);

#endif

