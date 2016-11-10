#ifndef SKA_ASTROACCELERATE_SAVE_DATA_H_
#define SKA_ASTROACCELERATE_SAVE_DATA_H_

//void save_data(cudaStream_t stream, float *device_pointer, float *host_pointer, size_t size);
void save_data(float *device_pointer, float *host_pointer, size_t size);

#endif

