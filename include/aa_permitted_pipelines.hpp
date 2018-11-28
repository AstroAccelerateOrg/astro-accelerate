//
//  aa_permitted_pipelines.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 22/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_PERMITTED_PIPELINES_HPP
#define ASTRO_ACCELERATE_PERMITTED_PIPELINES_HPP

#include "aa_compute.hpp"
#include "aa_permitted_pipelines.hpp"
#include "aa_permitted_pipelines_0.hpp"
#include "aa_permitted_pipelines_1.hpp"

namespace astroaccelerate {

  /**
   * This class is used to check whether a pipeline is valid and permitted.
   * This class can also be used to obtain a valid pipeline.
   */

  class aa_permitted_pipelines {
  public:
    
    //Example valid pipelines
    static const aa_compute::pipeline pipeline0;
    static const aa_compute::pipeline pipeline1;
    static const aa_compute::pipeline pipeline2;
    static const aa_compute::pipeline pipeline3;
    static const aa_compute::pipeline pipeline4;
    static const aa_compute::pipeline pipeline5;
    
    static bool is_permitted(const aa_compute::pipeline &pipeline) {
      /**
       * Check if pipeline exists and is valid.
       */
        
      if(pipeline == pipeline0) {
	return true;
      }
      else if(pipeline == pipeline1) {
	return true;
      }
      else if(pipeline == pipeline2) {
	return true;
      }
      else {
	return false;
      }
        
      return false;
    }
    
  };

} //namespace astroaccelerate

#endif /* ASTRO_ACCELERATE_PERMITTED_PIPELINES_HPP */
