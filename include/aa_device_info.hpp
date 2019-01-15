#ifndef ASTRO_ACCELERATE_AA_DEVICE_INFO_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_INFO_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

#include "aa_log.hpp"

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

    aa_device_info(const aa_device_info&) = delete; /** \brief Delete copy and move constructors. */
    aa_device_info(aa_device_info&&) = delete; /** \brief Delete copy and move constructors. */
    aa_device_info& operator=(const aa_device_info&) = delete; /** \brief Delete copy and move constructors. */
    aa_device_info& operator=(aa_device_info&&) = delete; /** \brief Delete copy and move constructors. */

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

    static aa_device_info& instance() {
      static aa_device_info m_instance;
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
	LOG(log_level::dev_debug, "cudaGetDeviceCount returned" + std::to_string(error_id) + "->" + cudaGetErrorString(error_id));
	return false;
      }
        
      // This function call returns 0 if there are no CUDA capable devices.
      if(deviceCount == 0) {
	LOG(log_level::notice, "There are no available device(s) that support CUDA");
      }
      else {
	LOG(log_level::notice, "Detected " +std::to_string(deviceCount) + " CUDA Capable device(s)");
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
	LOG(log_level::notice, "Device info " + std::to_string(i) + " free " + std::to_string(free));
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
	LOG(log_level::notice, "No card has yet been selected. Defaulting to card with ID 0 and proceeding.");
	selected_card_idx = 0;
        if(!check_for_devices()) {
          LOG(log_level::error, "Could not check devices.");
          return false;
        }
      }
    
      for(size_t i = 0; i < m_card_info.size(); i++) {
	if(m_card_info.at(i).card_number == id) {
#ifdef ASTRO_ACCELERATE_VERSION_H_DEFINED
	  std::string device_compute_capability = std::to_string(m_card_info.at(i).compute_capability_major) + std::to_string(m_card_info.at(i).compute_capability_minor);
	  if(std::to_string(ASTRO_ACCELERATE_CUDA_SM_VERSION) > device_compute_capability) {
	    LOG(log_level::error, "Compiled for compute capability " + std::to_string(ASTRO_ACCELERATE_CUDA_SM_VERSION) + ".\n" + "The requested device has capability " + device_compute_capability + ".");
	    return false;
	  }
	  else {
	    LOG(log_level::notice, "Application binary compiled for compute capability "+std::to_string(ASTRO_ACCELERATE_CUDA_SM_VERSION)+".");
	    LOG(log_level::notice, "The requested device has capability "+device_compute_capability+".");
	  }
#else
	  LOG(log_level::notice, "Because #include \"version.h\" is not created by this build system, the compute capability of the device cannot be determined.");
	  LOG(log_level::notice, "Please consider compiling using the CMakeLists file provided in the repository.");
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
	LOG(log_level::notice, "No card has yet been selected. Defaulting to card with ID 0 and proceeding.");
	selected_card_idx = 0;
	if(!check_for_devices()) {
	  LOG(log_level::error, "Could not check devices.");
	  return false;
	}
      }
      if(m_card_info.at(selected_card_idx).user_requested_memory_for_allocation + mem >= gpu_memory()) {
	LOG(log_level::error,
	    "Device reports that the additional requested memory ("+std::to_string(mem)+")"+", in addition to the already requested memory ("+std::to_string(requested())+"), would exceed the total free memory on the device ("+std::to_string((unsigned long long)gpu_memory())+").");
	return false;
      }
      m_card_info.at(selected_card_idx).user_requested_memory_for_allocation += mem;
      return true;
    }

    /** \returns The amount of memory requested on the currently selected, or 0 if cards have not been checked. */
    size_t requested() {
      if(!is_init) {
	LOG(log_level::notice, "No card has yet been selected. Defaulting to card with ID 0 and proceeding.");
	selected_card_idx = 0;
	if(!check_for_devices()) {
          LOG(log_level::error, "Could not check devices.");
          return false;
        }
	return 0;
      }
      return m_card_info.at(selected_card_idx).user_requested_memory_for_allocation;
    }

    /** \brief Reset all requested memory to 0 for a given card ID. */
    bool reset_requested_memory_on_card(const CARD_ID &id) {
      if(!is_init) {
	return false;
      }
      
      if((unsigned long)id < m_card_info.size()) {
	m_card_info.at(selected_card_idx).user_requested_memory_for_allocation = 0;
	return true;
      }
      
      return false;
    }

    /** \brief Static method for printing member data for an instance of aa_card_info. */
    static bool print_card_info(const aa_card_info &card) {
      LOG(log_level::dev_debug, "CARD INFORMATION:");
      LOG(log_level::dev_debug, "Name:\t\t\t"+card.name);
      LOG(log_level::dev_debug, "Compute capability:\t" + std::to_string(card.compute_capability_major) + "." + std::to_string(card.compute_capability_minor));
      LOG(log_level::dev_debug, "Clock rate:\t\t" + std::to_string(card.clock_rate));
      LOG(log_level::dev_debug, "Memory clock rate:\t" + std::to_string(card.memory_clock_rate));
      LOG(log_level::dev_debug, "Memory bus width:\t" + std::to_string(card.memory_bus_width));
      LOG(log_level::dev_debug, "ID # on this machine:\t" + std::to_string(card.card_number));
      LOG(log_level::dev_debug, "Driver version:\t\t" + std::to_string(card.driver_version));
      LOG(log_level::dev_debug, "l2_cache_size:\t\t" + std::to_string(card.l2_cache_size));
      LOG(log_level::dev_debug, "Global memory:\t\t" + std::to_string(card.global_memory));
      LOG(log_level::dev_debug, "Total memory:\t\t" + std::to_string(card.total_memory));
      LOG(log_level::dev_debug, "Free memory:\t\t" + std::to_string(card.free_memory));
      return true;
    }
    
  private:
    aa_device_info() : selected_card_idx(0) {

    }

     
    /** \brief The destructor is public so static references returned by the instance() method are deleted. */
    ~aa_device_info() {

    }
    
    static bool is_init; /**< Flag to indicate whether cards on the machine have been checked/initialised. */
    static std::vector<aa_card_info> m_card_info; /** Stores all card information for all cards on the machine. */
    CARD_ID selected_card_idx;  /**< Index into m_card_info for the current selected card. */
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_INFO_HPP
