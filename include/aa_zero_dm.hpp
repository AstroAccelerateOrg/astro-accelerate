//
//  aa_zero_dm.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 05/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_zero_dm_hpp
#define aa_zero_dm_hpp

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

#endif /* aa_zero_dm_hpp */
