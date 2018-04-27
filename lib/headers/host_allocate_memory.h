#ifndef ASTROACCELERATE_ALLOCATEMEMORY_H_
#define ASTROACCELERATE_ALLOCATEMEMORY_H_

void allocate_memory_cpu_input(unsigned short **input_buffer, size_t nsamp, size_t nchans);

void allocate_memory_cpu_output(float ****output_buffer, DDTR_Plan *DDTR_plan);

void allocate_memory_gpu(FILE **fp, size_t gpu_memory, int max_samps, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, unsigned short **input_buffer, float ****output_buffer, unsigned short **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize);

#endif
