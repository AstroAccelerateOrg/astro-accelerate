#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/sysinfo.h>
#include "aa_host_info.hpp"

namespace astroaccelerate {

void host_info(struct sysinfo *host_info)
{
        if ( sysinfo(host_info)== -1){
                printf("\n!!!Error on host system info!!!\n");
                return;
        }
}


// adding the possibility to read from a file /proc/meminfo; The sysinfo is showing free ram but without cached and buffers
// the line from meminfo with MemAvailable is the memory available to launch the application without touching the swap
int GetRamInKB(size_t *host_memory)
{
    FILE *meminfo = fopen("/proc/meminfo", "r");
    if(meminfo == NULL){
	printf("\n!!!Error on host system info!!!\n");
        return -1;
    }

    char line[256];
    while(fgets(line, sizeof(line), meminfo))
    {
       if(sscanf(line, "MemAvailable: %zu kB", host_memory) == 1)
        {
            fclose(meminfo);
//		printf("\n\t\t Ram: %d", *host_memory);
		*host_memory = *host_memory*1024;
            return *host_memory;
        } 
    }
	// return if not find the line
    fclose(meminfo);
    return -1;
}

void host_mem_error(unsigned int inputsize, unsigned int host_memory, const char *type)
{
	printf("\n\nCan't allocate %s memory of size: %u MiB. Host available memory only: %u MiB.\n",type, inputsize, host_memory);
	return;
}

} //namespace astroaccelerate
