//
//  aa_permitted_pipelines_1.hpp
//  aapipeline
//
//  Created by Cees Carels on Friday 02/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_PERMITTED_PIPELINES_1_HPP
#define ASTRO_ACCELERATE_PERMITTED_PIPELINES_1_HPP


#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <stdio.h>
#include "aa_ddtr_strategy.hpp"
#include "aa_ddtr_plan.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_device_load_data.hpp"
#include "aa_bin_gpu.hpp"
#include "aa_zero_dm.hpp"
#include "aa_zero_dm_outliers.hpp"
#include "aa_corner_turn.hpp"
#include "aa_dedisperse.hpp"

namespace astroaccelerate {

void run_pipeline_1(const aa_filterbank_metadata &metadata, const aa_ddtr_strategy &ddtr_strategy, unsigned short *const input_buffer);

}
  
#endif /* ASTRO_ACCELERATE_PERMITTED_PIPELINES_1_HPP */
