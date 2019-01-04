#ifndef ASTRO_ACCELERATE_AA_DEVICE_MEMORY_MANAGER_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_MEMORY_MANAGER_HPP

#include <stdio.h>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include "aa_device_info.hpp"

namespace astroaccelerate {

  /**
   * \class aa_device_memory_manager aa_device_memory_manager.hpp "include/aa_device_memory_manager.hpp"
   * \brief Class that tracks the amount of memory that modules have requested.
   * \details Modules can notify the manager that they will request memory, but this is not enforced.
   * \details This class exists to enable the user to configure plan and strategy objects without requiring that they allocate memory to hold large input/output data.
   * \details Therefore, it is the pipeline object that performs memory allocations for input/output data.
   * \warning Modules may choose not to use the manager, so it may not reflect the true amount of memory used, developers must implement this functionality consistently.
   * \warning Modules need to inform the manager when they de-allocate too.
   * \author Cees Carels.
   * \date 27 November 2018. 
   */

  class aa_device_memory_manager {
  public:

    /** \brief Constructor for aa_device_memory_manager. */
    aa_device_memory_manager() : m_init(false) {
      //Set the 
    }

    /** \brief Attempt to request memory. */
    bool request(const size_t &mem) {
      if(mem < 17179869184 && m_mem.size()) { // Manual cutoff point
	aa_device_info* device_info = aa_device_info::instance();
	m_mem.at(device_info->currently_selected_card_id()) += mem;
	return true;
      }
      else {
	return false;}
    }

    /** \returns The total memory planned to be allocated on the currently selected card. */
    size_t requested() {
      aa_device_info* device_info = aa_device_info::instance();
      return m_mem.size() ? m_mem.at(device_info->currently_selected_card_id()) : 0;
    }

    /** \returns The free memory on the currently selected card. */
    size_t free_memory() {
      aa_device_info* device_info = aa_device_info::instance();
      device_info->gpu_memory();
    }
  
  private:
    static bool m_init;
    static unsigned long m_mem;

  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_MEMORY_MANAGER_HPP
