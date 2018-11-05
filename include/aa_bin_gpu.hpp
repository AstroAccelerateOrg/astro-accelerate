//
//  aa_bin_gpu.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 05/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_bin_gpu_hpp
#define aa_bin_gpu_hpp

#include <stdio.h>
#include <vector_types.h>

#ifndef AA_CUDA
#define AA_CUDA 1
#endif

#include "params.hpp"
#include "gpu_timer.hpp"
#include "device_binning_kernel.hpp"
#include "device_corner_turn_kernel.hpp"

void bin_gpu(unsigned short *const d_input, float *const d_output, const int nchans, const int nsamp);
int GPU_DiT_v2_wrapper(float *d_input, float *d_output, int nDMs, int nTimesamples);

#endif /* aa_bin_gpu_hpp */
