//
//  aa_zero_dm_outliers.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 05/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_zero_dm_outliers_hpp
#define aa_zero_dm_outliers_hpp

#include <stdio.h>
#include <time.h>
#include <vector_types.h>

#include "device_zero_dm_outliers_kernel.hpp"

void zero_dm_outliers(unsigned short *const d_input, const int nchans, const int nsamp);

#endif /* aa_zero_dm_outliers_hpp */
