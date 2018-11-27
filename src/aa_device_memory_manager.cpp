#include "aa_device_memory_manager.hpp"

namespace astroaccelerate {
  // Declaration of static memory in source so that the reference can be resolved
  // when using class aa_device_memory_manager.hpp
  bool aa_device_memory_manager::m_init = false;
  std::vector<size_t> aa_device_memory_manager::m_mem;
}
