#ifndef __ALLOCATEMEMORY__
#define __ALLOCATEMEMORY__

void allocate_memory(FILE **fp, size_t gpu_memory, int max_samps, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, float **input_buffer, float ****output_buffer, float **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize);

#endif
