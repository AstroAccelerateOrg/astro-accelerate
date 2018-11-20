//
//  aa_analysis_plan.hpp
//  aapipeline
//
//  Created by Cees Carels on Tuesday 23/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_ANALYSIS_PLAN_HPP
#define ASTRO_ACCELERATE_ANALYSIS_PLAN_HPP

#include <stdio.h>

#include "aa_ddtr_plan.hpp"

namespace astroaccelerate {

  class aa_analysis_plan {
  public:

    enum class selectable_candidate_algorithm : int {
	off = 0,
	on			  
    };
    
    aa_analysis_plan() {
      
    }
    
    aa_analysis_plan(const aa_ddtr_strategy &ddtr_strategy,
		     const float &sigma_cutoff,
		     const float &sigma_constant,
		     const float &max_boxcar_width_in_sec,
		     const selectable_candidate_algorithm &candidate_algorithm,
		     const bool &enable_sps_baseline_noise) : m_ddtr_strategy(ddtr_strategy),
							      m_sigma_cutoff(sigma_cutoff),
							      m_sigma_constant(sigma_constant),
							      m_max_boxcar_width_in_sec(max_boxcar_width_in_sec),
							      m_candidate_algorithm(candidate_algorithm),
							      m_enable_sps_baseline_noise(enable_sps_baseline_noise) {
      
    }
    
    const aa_ddtr_strategy ddtr_strategy() const {
      return m_ddtr_strategy;
    }

    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }

    float sigma_constant() const {
      return m_sigma_constant;
    }

    float max_boxcar_width_in_sec() const {
      return m_max_boxcar_width_in_sec;
    }

    aa_analysis_plan::selectable_candidate_algorithm candidate_algorithm() const {
      return m_candidate_algorithm;
    }

    bool enable_sps_baseline_noise() const {
      return m_enable_sps_baseline_noise;
    }

  private:
    aa_ddtr_strategy    m_ddtr_strategy;
    float               m_sigma_cutoff;
    float               m_sigma_constant;
    float               m_max_boxcar_width_in_sec;
    aa_analysis_plan::selectable_candidate_algorithm m_candidate_algorithm;
    bool                m_enable_sps_baseline_noise;
  };
  
} //namespace astroaccelerate

#endif /* ASTRO_ACCELERATE_ANALYSIS_PLAN_HPP */
