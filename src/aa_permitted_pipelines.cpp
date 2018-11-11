//
//  aa_permitted_pipelines.cpp
//  aapipeline
//
//  Created by Cees Carels on Monday 22/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_permitted_pipelines.hpp"

namespace astroaccelerate {
  
const aa_compute::pipeline aa_permitted_pipelines::pipeline0 = {aa_compute::modules::empty};
const aa_compute::pipeline aa_permitted_pipelines::pipeline1 = {aa_compute::modules::dedispersion};
const aa_compute::pipeline aa_permitted_pipelines::pipeline2 = {aa_compute::modules::dedispersion, aa_compute::modules::analysis};

} //namespace astroaccelerate
