//
//  aa_device_load_data.cpp
//  aapipeline
//
//  Created by Cees Carels on Friday 02/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_device_load_data.hpp"
#include "device_dedispersion_kernel.hpp"
#include "device_SPS_inplace_kernel.hpp"

#include <helper_cuda.h>
#include <iostream>

void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short *host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts) {
    if(i == -1) {
        long int length = ( t_processed + maxshift );
        size_t size = nchans * length * sizeof(unsigned short);
        cudaMemcpy(device_pointer, host_pointer, size, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());
        set_device_constants_dedispersion_kernel(nchans, length, t_processed, dmshifts);
	checkCudaErrors(cudaGetLastError());
    }
    else if(i > 0) {
        long int length = ( t_processed + maxshift );
        set_device_constants_dedispersion_kernel(length, t_processed);
    }

    float h_sqrt_taps[PD_MAXTAPS + 1];
    for(int f = 0; f <= PD_MAXTAPS; f++) {
      h_sqrt_taps[f] = (float) sqrt((double) f);
    }
    std::cout << "Checking for sqrt_taps" << std::endl;
    cudaError_t myError = cudaMemcpyToSymbol(c_sqrt_taps, &h_sqrt_taps, ( PD_MAXTAPS + 1 ) * sizeof(float));
    if(myError == cudaSuccess) {
      std::cout << "Success " << std::endl;
    }
    else {
      std::cout << "Not success" << std::endl;
    }
    std::cout << "The error was " << cudaGetErrorName (myError) << std::endl;

}
