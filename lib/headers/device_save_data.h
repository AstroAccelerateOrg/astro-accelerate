#ifndef ASTROACCELERATE_SAVE_DATA_H_
#define ASTROACCELERATE_SAVE_DATA_H_

//void save_data(cudaStream_t stream, float *device_pointer, float *host_pointer, size_t size);
void save_data(float *device_pointer, float *host_pointer, size_t size);
void save_data_offset(float *host_pointer_pinned, int device_offset, float *host_pointer, int host_offset, size_t size, cudaStream_t streams);

#endif

