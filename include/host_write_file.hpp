#ifndef ASTRO_ACCELERATE_HOST_WRITE_FILE_HPP
#define ASTRO_ACCELERATE_HOST_WRITE_FILE_HPP

void write_output(int i, int t_processed, int ndms, size_t  gpu_memory, float *output_buffer, size_t gpu_outputsize, float *dm_low, float *dm_high);

#endif
