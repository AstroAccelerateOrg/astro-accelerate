#ifndef ASTRO_ACCELERATE_AA_DEVICE_MEMORY_MANAGER_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_MEMORY_MANAGER_HPP

/**
 * Class that tracks the total amount of memory that various modules have requested.
 */

#include <stdio.h>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

namespace astroaccelerate {

class aa_device_memory_manager {
public:
  bool request(const size_t &mem) {
    if(!m_init) {  
      if(!init()) {
	return false;
      }
    }
    
    int device_id = 0;
    cudaError_t error_id = cudaGetDevice(&device_id);
    if(error_id != cudaSuccess) {
      printf("cudaGetDevice returned %d\n-> %s\n", (int) error_id, cudaGetErrorString(error_id));
      printf("Result = FAIL\n");
      return false;
    }
    
    if(mem > 17179869184) { // Manual cutoff point
      return false;
    }
    
    if(m_mem.size()) {
      m_mem.at(device_id) += mem;
    }
    return true;
  }

  size_t requested() {
    /**
     * Total memory planned to be allocated.
     */
    if(!m_init) {
      if(!init()) {
        return false;
      }
    }
    
    int	device_id = 0;
    cudaError_t error_id = cudaGetDevice(&device_id);
    if(error_id != cudaSuccess) {
      printf("cudaGetDevice returned %d\n-> %s\n", (int) error_id, cudaGetErrorString(error_id));
      printf("Result = FAIL\n");
      return false;
    }
    
    return m_mem.size() ? m_mem.at(device_id) : 0;
  }

  size_t free_memory() {
    /**
     * Free memory on the currently selected card.
     * TODO: Make a distinction between memory that is requested but not yet
     * allocated, and the actual amount of memory that is free.
     */
    if(!m_init) {
      if(!init()) {
        return false;
      }
    }
    
    int device_id = 0;
    cudaError_t error_id = cudaGetDevice(&device_id);
    if(error_id != cudaSuccess) {
      printf("cudaGetDevice returned %d\n-> %s\n", (int) error_id, cudaGetErrorString(error_id));
      printf("Result = FAIL\n");
      return false;
    }
    
    size_t free = 0;
    size_t total = 0;
    error_id = cudaMemGetInfo(&free, &total);
    if(error_id != cudaSuccess) {
      printf("cudaGetMemInfo returned %d\n-> %s\n", (int) error_id, cudaGetErrorString(error_id));
      printf("Result = FAIL\n");
    }
    
    return (m_mem.at(device_id) > free) ? free : free-m_mem.at(device_id);
  }
  
private:
  static bool m_init;
  static std::vector<size_t> m_mem;

  bool init() {
    int deviceCount = 0;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
    
    if(error_id != cudaSuccess) {
      printf("cudaGetDeviceCount returned %d\n-> %s\n", (int) error_id, cudaGetErrorString(error_id));
      printf("Result = FAIL\n");
      return false;
    }
    
    // This function call returns 0 if there are no CUDA capable devices. 
    if(deviceCount == 0) {
      printf("There are no available device(s) that support CUDA\n");
      return false;
    }
    else {
      printf("Detected %d CUDA Capable device(s)\n", deviceCount);
      m_mem.resize(deviceCount, 0);
    }
    m_init = true;
    return true;
  }
  
};

} //namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_MEMORY_MANAGER_HPP
