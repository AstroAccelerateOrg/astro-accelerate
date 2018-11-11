//
//  aa_zero_dm_outliers.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 05/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_ZERO_DM_OUTLIERS_HPP
#define ASTRO_ACCELERATE_ZERO_DM_OUTLIERS_HPP

#include <stdio.h>
#include <time.h>
#include <vector_types.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_zero_dm_outliers_kernel.hpp"

namespace astroaccelerate {

void zero_dm_outliers(unsigned short *const d_input, const int nchans, const int nsamp);

} //namespace astroaccelerate
  
#endif /* ASTRO_ACCELERATE_ZERO_DM_OUTLIERS_HPP */
