#ifndef ASTRO_ACCELERATE_AA_DEVICE_INFO_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_INFO_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

/**
 * If the build system (CMake) defines aa_version.hpp, then include it.
 */
#ifdef ASTRO_ACCELERATE_VERSION_H_DEFINED
#include "aa_version.hpp"
#endif

namespace astroaccelerate {

  /**
   * \class aa_device_info aa_device_info.hpp "include/aa_device_info.hpp" 
   * \brief Obtain information about available GPUs and select the GPU to use for data processing.
   * \author Cees Carels.
   * \date 1 November 2018.
   */

  class aa_device_info {
  public:
    typedef int CARD_ID;

    /**
     * \struct aa_card_info
     * \brief Struct to contain CUDA card information.
     */
    struct aa_card_info {
      int compute_capability_major;
      int compute_capability_minor;
      float clock_rate;
      float memory_clock_rate;
      float memory_bus_width;
      CARD_ID card_number;
      int driver_version;
      int runtime_version;
      int l2_cache_size;
      size_t global_memory;
      size_t total_memory;
      size_t free_memory;
      size_t user_requested_memory_for_allocation; // Memory the user wishes to allocate, may or may not already have been allocated.
      std::string name;
    };

    static aa_device_info* instance() {
      if(!m_instance) {
	m_instance = new aa_device_info;
      }
      return m_instance;
    }

    /** \brief Checks for GPUs on the machine.
     * \details This is an overloaded function without parameters.
     * \details All card information is contained inside the member variable m_card_info.
     * \returns A boolean to indicate whether this operation was successful.
     */
    bool check_for_devices() {
      return check_for_devices(m_card_info);
    }
    
    /** \brief Checks for GPUs on the machine.
     * \details Provides the information back to the user via the card_info parameters.
     * \returns A boolean to indicate whether this operation was successful.
     */
    bool check_for_devices(std::vector<aa_card_info> &card_info) {
      m_card_info.clear(); //Clear the vector, so that in case check_for_devices is called multiple times, it does not push_back multiple entries
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
      }
      else {
	printf("Detected %d CUDA Capable device(s)\n", deviceCount);
      }
        
      for(int i = 0; i < deviceCount; i++) {
	aa_card_info tmp;
	tmp.card_number = i;
            
	cudaSetDevice(i);
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, i);
            
	int driverVersion = 0;
	int runtimeVersion = 0;
	cudaDriverGetVersion(&driverVersion);
	cudaRuntimeGetVersion(&runtimeVersion);
            
	tmp.name = deviceProp.name;
	tmp.driver_version = driverVersion;
	tmp.runtime_version = runtimeVersion;
	tmp.global_memory = deviceProp.totalGlobalMem;
	tmp.clock_rate = deviceProp.clockRate;
	tmp.memory_clock_rate = deviceProp.memoryClockRate;
	tmp.memory_bus_width = deviceProp.memoryBusWidth;
	tmp.l2_cache_size = deviceProp.l2CacheSize;
	tmp.compute_capability_major = deviceProp.major;
	tmp.compute_capability_minor = deviceProp.minor;
	//All additional properties read out here...
            
	size_t free     = 0;
	size_t total    = 0;
	cudaMemGetInfo(&free, &total);
	tmp.free_memory = free;
	tmp.total_memory = total;
	std::cout << "Device info " << i << " free " << free << std::endl;
	tmp.user_requested_memory_for_allocation = 0;
	m_card_info.push_back(std::move(tmp));
      }
        
      selected_card_idx = 0;   // The default selected card is the first one.
      card_info = m_card_info; // Assign the m_card_info to the parameter passed.
      is_init = true;
      return true;
    }
    
    /** \returns A boolean to indicate whether selecting the card was successful. */
    bool init_card(const CARD_ID &id, aa_card_info &card_info) {
      if(!is_init) {
	std::cout << "NOTICE: No card has yet been selected. Defaulting to card with ID 0 and proceeding." << std::endl;
        if(!check_for_devices()) {
          std::cout << "ERROR:  Could not check devices." << std::endl;
          return false;
        }
      }
    
      for(size_t i = 0; i < m_card_info.size(); i++) {
	if(m_card_info.at(i).card_number == id) {
#ifdef ASTRO_ACCELERATE_VERSION_H_DEFINED
	  std::string device_compute_capability = std::to_string(m_card_info.at(i).compute_capability_major) + std::to_string(m_card_info.at(i).compute_capability_minor);
	  if(ASTRO_ACCELERATE_CUDA_SM_VERSION > device_compute_capability) {
	    std::cout << "ERROR: Compiled for compute capability " << ASTRO_ACCELERATE_CUDA_SM_VERSION << "." << std::endl;
	    std::cout << "The requested device has capability " << device_compute_capability << "." << std::endl;
	    return false;
	  }
	  else {
	    std::cout << "NOTICE: Application binary compiled for compute capability " << ASTRO_ACCELERATE_CUDA_SM_VERSION << "." << std::endl;
	    std::cout << "        The requested device has capability " << device_compute_capability << "." << std::endl;
	  }
#else
	  std::cout << "NOTICE: Because #include \"version.h\" is not created by this build system, the compute capability of the device cannot be determined." << std::endl;
	  std::cout << "        Please consider compiling using the CMakeLists file provided in the repository." << std::endl;
#endif
	  selected_card_idx = i;
	  cudaSetDevice(i);
	  size_t free     = 0;
	  size_t total    = 0;
	  cudaMemGetInfo(&free, &total);
	  m_card_info.at(i).free_memory = free;
	  m_card_info.at(i).total_memory = total;
	  card_info = m_card_info.at(i);
	  return true;
	}
      }
        
      return true;
    }
    
    /** \returns The currently available free memory on the currently selected GPU as reported by the CUDA driver. */
    size_t gpu_memory() const {
      if(is_init) {
	return m_card_info.at(selected_card_idx).free_memory;
      }
        
      return 0;
    }

    CARD_ID currently_selected_card_id() const {
      return selected_card_idx;
    }

    /**
     * \brief Request memory on the currently selected card.
     * \returns True if the request was successful, false otherwise (e.g. cards have not been checked).
     */
    bool request(const size_t mem) {
      if(!is_init) {
	std::cout << "NOTICE: No card has yet been selected. Defaulting to card with ID 0 and proceeding." << std::endl;
	if(!check_for_devices()) {
	  std::cout << "ERROR:  Could not check devices." << std::endl; 
	  return false;
	}
      }
      if(m_card_info.at(selected_card_idx).user_requested_memory_for_allocation + mem >= gpu_memory()) {
	std::cout << "ERROR:  Device reports that the additional requested memory (" << mem << "), "
		  << "in addition to the already requested memory (" << requested() << "), would exceed the total free memory on the device (" << (unsigned long long)gpu_memory() << ")."
		  << std::endl;
	return false;
      }
      m_card_info.at(selected_card_idx).user_requested_memory_for_allocation += mem;
      return true;
    }

    /** \returns The amount of memory requested on the currently selected, or 0 if cards have not been checked. */
    size_t requested() {
      if(!is_init) {
	std::cout << "NOTICE: No card has yet been selected. Defaulting to card with ID 0 and proceeding." << std::endl;
	if(!check_for_devices()) {
          std::cout << "ERROR:  Could not check devices." << std::endl;
          return false;
        }
	return 0;
      }
      return m_card_info.at(selected_card_idx).user_requested_memory_for_allocation;
    }

    /** \brief Static method for printing member data for an instance of aa_card_info. */
    static bool print_card_info(const aa_card_info &card) {
      std::cout << "CARD INFORMATION:" << std::endl;
      std::cout << "Name:\t\t\t" << card.name << std::endl;
      std::cout << "Compute capability:\t"
		<< card.compute_capability_major << "."
		<< card.compute_capability_minor << std::endl;
      std::cout << "Clock rate:\t\t" << card.clock_rate << std::endl;
      std::cout << "Memory clock rate:\t" << card.memory_clock_rate << std::endl;
      std::cout << "Memory bus width:\t" << card.memory_bus_width << std::endl;
      std::cout << "ID # on this machine:\t" << card.card_number << std::endl;
      std::cout << "Driver version:\t\t" << card.driver_version << std::endl;
      std::cout << "l2_cache_size:\t\t" << card.l2_cache_size << std::endl;
      std::cout << "Global memory:\t\t" << card.global_memory << std::endl;
      std::cout << "Total memory:\t\t" << card.total_memory << std::endl;
      std::cout << "Free memory:\t\t" << card.free_memory << std::endl;
      return true;
    }
  
  private:
    /** \brief The only constructor is private so that only the singleton instance pointer can be used to access the class. */
    aa_device_info() : selected_card_idx(0) {
      
    }
    
    static aa_device_info* m_instance;
    static bool is_init; /**< Flag to indicate whether cards on the machine have been checked/initialised. */
    static std::vector<aa_card_info> m_card_info; /** Stores all card information for all cards on the machine. */
    CARD_ID selected_card_idx;  /**< Index into m_card_info for the current selected card. */
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_INFO_HPP
