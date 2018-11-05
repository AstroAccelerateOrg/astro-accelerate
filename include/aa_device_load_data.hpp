//
//  aa_device_load_data.hpp
//  aapipeline
//
//  Created by Cees Carels on Friday 02/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_device_load_data_hpp
#define aa_device_load_data_hpp

#include <stdio.h>
#include <math.h>

void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short *host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts);

#endif /* aa_device_load_data_hpp */
