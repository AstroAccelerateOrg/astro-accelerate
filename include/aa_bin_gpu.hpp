//
//  aa_bin_gpu.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 05/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_BIN_GPU_HPP
#define ASTRO_ACCELERATE_BIN_GPU_HPP

#include <stdio.h>
#include <vector_types.h>

#include "params.hpp"
#include "gpu_timer.hpp"
#include "device_binning_kernel.hpp"
#include "device_corner_turn_kernel.hpp"

void bin_gpu(unsigned short *const d_input, float *const d_output, const int nchans, const int nsamp);
int GPU_DiT_v2_wrapper(float *d_input, float *d_output, int nDMs, int nTimesamples);

#endif /* ASTRO_ACCELERATE_BIN_GPU_HPP */
