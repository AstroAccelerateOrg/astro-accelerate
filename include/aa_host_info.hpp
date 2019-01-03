#ifndef ASTRO_ACCELERATE_AA_HOST_INFO_HPP
#define ASTRO_ACCELERATE_AA_HOST_INFO_HPP

#include <sys/sysinfo.h>

namespace astroaccelerate {
  /** \brief Function that gets sysinfo host information. */
  void host_info(struct sysinfo *host_info);

  /** \brief Function that prints memory errors. */
  void host_mem_error(unsigned int inputsize, unsigned int host_memory, const char *type);

  /** \brief Function that prints host memory in KB. */
  int GetRamInKB(size_t *host_memory);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_HOST_INFO_HPP
