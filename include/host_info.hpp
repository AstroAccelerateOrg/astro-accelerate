#ifndef ASTRO_ACCELERATE_HOST_INFO_HPP
#define ASTRO_ACCELERATE_HOST_INFO_HPP

#include <sys/sysinfo.h>

namespace astroaccelerate {

void host_info(struct sysinfo *host_info);
void host_mem_error(unsigned int inputsize, unsigned int host_memory, const char *type);
int GetRamInKB(size_t *host_memory);

} //namespace astroaccelerate
  
#endif
