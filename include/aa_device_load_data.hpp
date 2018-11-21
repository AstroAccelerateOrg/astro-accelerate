//
//  aa_device_load_data.hpp
//  aapipeline
//
//  Created by Cees Carels on Friday 02/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_DEVICE_LOAD_DATA_HPP
#define ASTRO_ACCELERATE_DEVICE_LOAD_DATA_HPP

#include <stdio.h>
#include <math.h>

namespace astroaccelerate {

void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short const*const host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts);

} //namespace astroaccelerate

#endif /* ASTRO_ACCELERATE_DEVICE_LOAD_DATA_HPP */
