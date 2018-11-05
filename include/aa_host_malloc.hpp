//
//  aa_host_malloc.hpp
//  aapipeline
//
//  Created by Cees Carels on Wednesday 31/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_host_malloc_hpp
#define aa_host_malloc_hpp

#include <stdio.h>

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

class aa_malloc {
public:
    template <typename T>
    bool aa_mem_malloc(T *ptr, const size_t &size) {
        if(size == 0) {
            //Cannot allocate memory of size 0.
            return false;
        }
        
        size_t host_memory = 0;
        if(GetRamInKB(host_memory) == -1) {
            return false;
        }
        
        if(host_memory < size) {
            return false;
        }
        else {
            ptr = (T*)malloc(size);
        }
        
        return true;
    }
    
    template <typename T>
    bool aa_mem_malloc(T **ptr, const size_t &size) {
        if(size == 0) {
            //Cannot allocate memory of size 0.
            return false;
        }
        
        size_t host_memory = 0;
        if(GetRamInKB(host_memory) == -1) {
            return false;
        }
        
        if(host_memory < size) {
            return false;
        }
        else {
            *ptr = (T*)malloc(size);
        }
        
        return true;
    }
    
    const size_t host_memory() const {
        size_t host_memory = 0;
        
        if(GetRamInKB(host_memory) == -1) {
            return false;
        }
        
        return host_memory;
    }
    
protected:
    
    const int GetRamInKB(size_t &host_memory) const {
#ifdef __APPLE__
        int mib[2];
        int64_t physical_memory;
        size_t length;
        
        // Get the Physical memory size
        mib[0] = CTL_HW;
        mib[1] = HW_MEMSIZE;
        length = sizeof(int64_t);
        sysctl(mib, 2, &physical_memory, &length, NULL, 0);
        host_memory = (size_t)(physical_memory);
        return 0;
#endif
        
#ifdef __linux__
        FILE *meminfo = fopen("/proc/meminfo", "r");
        if(meminfo == NULL){
            printf("ERROR: Could not access /proc/meminfo to read host memory.\n");
            return -1;
        }
        
        char line[256];
        while(fgets(line, sizeof(line), meminfo)) {
            if(sscanf(line, "MemAvailable: %zu kB", &host_memory) == 1) {
                fclose(meminfo);
                //        printf("\n\t\t Ram: %d", *host_memory);
                host_memory *= 1024;
                return 0;
            }
        }
        // return if not find the line
        fclose(meminfo);
        return -1;
#endif
        return -1;
    }
};

#endif /* aa_host_malloc_hpp */
