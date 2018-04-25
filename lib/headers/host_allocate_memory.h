#ifndef ASTROACCELERATE_ALLOCATEMEMORY_H_
#define ASTROACCELERATE_ALLOCATEMEMORY_H_

void allocate_memory_cpu_input(FILE **fp, size_t gpu_memory, int max_samps, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, unsigned short **input_buffer, float ****output_buffer, unsigned short **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize);

void allocate_memory_cpu_output(FILE **fp, size_t gpu_memory, int max_samps, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, unsigned short **input_buffer, float ****output_buffer, unsigned short **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize);

void allocate_memory_gpu(FILE **fp, size_t gpu_memory, int max_samps, int num_tchunks, int max_ndms, int total_ndms, int nsamp, int nchans, int nbits, int range, int *ndms, int **t_processed, unsigned short **input_buffer, unsigned short **input_buffer_small, float ****output_buffer, float **output_buffer_small, unsigned short **d_input, float **d_output, size_t *gpu_inputsize, size_t *gpu_outputsize, size_t *inputsize, size_t *outputsize);

void allocate_memory_MSD(float **d_MSD_workarea, unsigned short **d_MSD_output_taps, float **d_MSD_interpolated, float **h_MSD_interpolated, float **h_MSD_DIT, int **gmem_peak, int **temp_peak, unsigned long int MSD_maxtimesamples, int h_MSD_DIT_width, int nTimeSamples, size_t MSD_profile_size);

#endif
