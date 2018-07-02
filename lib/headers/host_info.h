#include <sys/sysinfo.h>

void host_info(struct sysinfo *host_info);
void host_mem_error(unsigned int inputsize, unsigned int host_memory, const char *type);
void GetRamInKB(size_t *host_memory);
