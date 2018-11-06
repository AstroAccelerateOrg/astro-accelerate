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

#include "params.hpp"

#include <helper_cuda.h>
#include <iostream>

void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short *host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts) {
    if(i == -1) {
        long int length = ( t_processed + maxshift );
        size_t size = (size_t)nchans * (size_t)length * (size_t)sizeof(unsigned short);
        cudaMemcpy(device_pointer, host_pointer, size, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());
        set_device_constants_dedispersion_kernel(nchans, length, t_processed, dmshifts);
    }
    else if(i > 0) {
        long int length = ( t_processed + maxshift );
        set_device_constants_dedispersion_kernel(length, t_processed);
    }

    float h_sqrt_taps[PD_MAXTAPS + 1];
    for(int f = 0; f <= PD_MAXTAPS; f++) {
      h_sqrt_taps[f] = (float) sqrt((double) f);
    }
    cudaMemcpyToSymbol(c_sqrt_taps, &h_sqrt_taps, ( PD_MAXTAPS + 1 ) * sizeof(float));
}
