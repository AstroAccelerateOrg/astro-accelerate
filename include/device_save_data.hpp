#ifndef ASTRO_ACCELERATE_DEVICE_SAVE_DATA_HPP
#define ASTRO_ACCELERATE_DEVICE_SAVE_DATA_HPP

// void save_data(cudaStream_t stream, float *device_pointer, float
// *host_pointer, size_t size);
void save_data(float *device_pointer, float *host_pointer, size_t size);
void save_data_offset(float *device_pointer,
                      int    device_offset,
                      float *host_pointer,
                      int    host_offset,
                      size_t size);

#endif
