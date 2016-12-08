#ifndef SKA_ASTROACCELERATE_ALLOCATEMEMORY_H_
#define SKA_ASTROACCELERATE_ALLOCATEMEMORY_H_

void allocate_memory_cpu_input(int nsamp, int nchans, unsigned short **input_buffer, size_t *inputsize);
void allocate_memory_cpu_output(int num_tchunks, int range, int *ndms, int **t_processed,
                                float ****output_buffer, float **d_output, size_t *outputsize);
void allocate_memory_gpu(int max_samps, int max_ndms, int nchans, int **t_processed,
                         unsigned short **d_input, float **d_output, size_t *gpu_inputsize,
                         size_t *gpu_outputsize);

#endif
