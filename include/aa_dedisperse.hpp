//
//  aa_dedisperse.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 05/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_dedisperse_hpp
#define aa_dedisperse_hpp

#include <stdio.h>
#include <math.h>
#include <vector_types.h>

#include <cuda.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>

#include "params.hpp"
#include "device_dedispersion_kernel.hpp"

#ifndef AA_CUDA
#define AA_CUDA 1
#endif

void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, unsigned short *d_input, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int const*const ndms, int nbits, int failsafe);

#endif /* aa_dedisperse_hpp */
