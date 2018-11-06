//
//  aa_zero_dm.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 05/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_ZERO_DM_HPP
#define ASTRO_ACCELERATE_ZERO_DM_HPP

#include <time.h>
#include <math.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>
#include <vector_types.h>

#include "params.hpp"
#include "device_zero_dm_kernel.hpp"

void zero_dm(unsigned short *const d_input, const int nchans, const int nsamp, const int nbits);

#endif /* ASTRO_ACCELERATE_ZERO_DM_HPP */
