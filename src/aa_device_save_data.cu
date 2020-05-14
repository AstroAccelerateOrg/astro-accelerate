#include <iostream>
#include <stdio.h>

namespace astroaccelerate {
  /** \brief Copy data and set up the GPU constants/variables. */
  void save_data(float *device_pointer, float *host_pointer, size_t size) {
    cudaMemcpy(host_pointer, device_pointer, size, cudaMemcpyDeviceToHost);
  }

  /** \brief Copy data and set up the GPU constants/variables. */
  void save_data_offset(float *device_pointer, size_t device_offset, float *host_pointer, size_t host_offset, size_t size) {
	cudaMemcpy(host_pointer + host_offset, device_pointer + device_offset, size, cudaMemcpyDeviceToHost);
  }
} //namespace astroaccelerate
