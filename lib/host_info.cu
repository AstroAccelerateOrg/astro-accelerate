#include <stdio.h>
#include <string.h>
#include <sys/sysinfo.h>

void host_info(struct sysinfo *host_info)
{
        if ( sysinfo(host_info)== -1){
                printf("\n!!!Error on host system info!!!\n");
                exit(0);
        }

//        printf("\nMemory: %zu MiB;\nMemory units: %zu;", (size_t)info.totalram/1024/1024, (size_t)info.mem_unit);
//        printf("\nAvailable memory: %zu MiB\n", (size_t)info.freeram/1024/1024);

}

void host_mem_error(unsigned int inputsize, unsigned int host_memory, const char *type)
{
	printf("\n\nCan't allocate %s memory of size: %u MiB. Host available memory only: %u MiB.\n",type, inputsize, host_memory);
}


//int main(int argc, char *argv[])
//{	
//	struct sysinfo info;
//	host_info(&info);
//	printf("\nAvailable memory: %zu MiB\n", (size_t)info.freeram/1024/1024);
//	
//	return 0;
//}


