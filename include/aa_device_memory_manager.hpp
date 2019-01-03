#ifndef ASTRO_ACCELERATE_AA_DEVICE_MEMORY_MANAGER_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_MEMORY_MANAGER_HPP

#include <stdio.h>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

namespace astroaccelerate {

  /**
   * \class aa_device_memory_manager aa_device_memory_manager.hpp "include/aa_device_memory_manager.hpp"
   * \brief Class that tracks the amount of memory that modules have requested.
   * \details Modules can notify the manager that they will request memory, but this is not enforced.
   * \details This class exists to enable the user to configure plan and strategy objects without requiring that they allocate memory to hold large input/output data.
   * \details Therefore, it is the pipeline object that performs memory allocations for input/output data.
   * \warning Modules may choose not to use the manager, so it may not reflect the true amount of memory used, developers must implement this functionality consistently.
   * \warning Modules need to inform the manager when they de-allocate too.
   * \todo The device ID must be made a parameter of the request method, so that the currently selected card is used, rather than defaulting to the 0th card in the machine.
   * \todo The device_info class may need a static member variable so that the currently selected GPU card ID can always be obtained from a static instance.
   * \author Cees Carels.
   * \date 27 November 2018. 
   */

  class aa_device_memory_manager {
  public:

    /**
     * Attempt to request memory.
     * \todo The GPU card ID must be specified, since at the moment it always defaults to allocating on card ID 0.
     */
    bool request(const size_t &mem) {
      if(!m_init) {  
	if(!init()) {
	  /** If there are no cards, then no memory can be requested. */
	  return false;
	}
      }
    
      int device_id = 0; /** \todo This must be made a parameter of the request method */
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
       * \todo Add a parameter to select the card ID.
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

    /**
     * Method that checks whether there are available GPUs.
     */
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

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_MEMORY_MANAGER_HPP
