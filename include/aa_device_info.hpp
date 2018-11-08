//
//  aa_device_init.hpp
//  aapipeline
//
//  Created by Cees Carels on Thursday 01/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_DEVICE_INFO_HPP
#define ASTRO_ACCELERATE_DEVICE_INFO_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>
#include <vector>
#include <string>

class aa_device_info {
public:
    typedef int CARD_ID;
    struct aa_card_info {
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
        std::string name;
    };
    
    bool check_for_devices() {
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
            //All additional properties read out here...
            
            size_t free     = 0;
            size_t total    = 0;
            cudaMemGetInfo(&free, &total);
	    tmp.free_memory = free;
	    tmp.total_memory = total;
	    std::cout << "Device info " << i << " free " << free << std::endl;
            m_card_info.push_back(std::move(tmp));
        }
        
        selected_card_idx = 0; // The default selected card is the first one.
        
        is_init = true;
        return true;
    }
    
    bool init_card(const CARD_ID &id, aa_card_info &card_info) {
        if(!is_init) {
            return false;
        }
        
        for(size_t i = 0; i < m_card_info.size(); i++) {
            if(m_card_info.at(i).card_number == id) {
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
    
    size_t gpu_memory() const {
        if(is_init) {
            return m_card_info.at(selected_card_idx).free_memory;
        }
        
        return 0;
    }
    
private:
    bool is_init;
    std::vector<aa_card_info> m_card_info;
    size_t selected_card_idx;  //Index into card_info for the current selected card.
};

#endif /* ASTRO_ACCELERATE_DEVICE_INFO_HPP */
