//
//  aa_device_corner_turn.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 05/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_device_corner_turn_hpp
#define aa_device_corner_turn_hpp

#include <stdio.h>
#include <time.h>
#include <vector_types.h>

#include "params.hpp"
#include "device_corner_turn_kernel.hpp"

#include <cuda.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>

#ifndef AA_CUDA
#define AA_CUDA 1
#endif

void corner_turn(unsigned short *const d_input, float *const d_output, const int nchans, const int nsamp);
int corner_turn(float *const d_input, float *const d_output, const int primary_size, const int secondary_size);
int corner_turn_SM(float *const d_input, float *const d_output, const int primary_size, const int secondary_size);

#endif /* aa_device_corner_turn_hpp */
