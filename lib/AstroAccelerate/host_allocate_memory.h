#ifndef SKA_ASTROACCELERATE_ALLOCATEMEMORY_H_
#define SKA_ASTROACCELERATE_ALLOCATEMEMORY_H_

#include "sps/InputData.h"
#include "sps/DedispersionPlan.h"

void allocate_memory_cpu_input(int nsamp, int nchans, unsigned short **input_buffer, size_t *inputsize);
void allocate_memory_cpu_input(ska::astroaccelerate::sps::DedispersionPlan const&, ska::astroaccelerate::sps::InputData const&);
void allocate_memory_cpu_output(size_t gpu_memory, int max_samps, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, unsigned short **input_buffer, float ****output_buffer, unsigned short **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize);

void allocate_memory_gpu(size_t gpu_memory, int max_samps, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, unsigned short **input_buffer, float ****output_buffer, unsigned short **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize);

#endif
