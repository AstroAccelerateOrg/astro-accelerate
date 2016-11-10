#ifndef SKA_ASTROACCELERATE_WRITE_H_
#define SKA_ASTROACCELERATE_WRITE_H_

void write_output(int i, int t_processed, int ndms, size_t  gpu_memory, float *output_buffer, size_t gpu_outputsize, float *dm_low, float *dm_high);

#endif
