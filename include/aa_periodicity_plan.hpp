//
//  aa_periodicity_plan.hpp
//  aapipeline
//
//  Created by Cees Carels on Tuesday 23/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP
#define ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP

#include <stdio.h>

namespace astroaccelerate {

  class aa_periodicity_plan {
  public:
    aa_periodicity_plan() {
      
    }
    
    aa_periodicity_plan(const float &sigma_cutoff,
			const float &OR_sigma_multiplier,
			const int   &nHarmonics,
			const int   &export_powers,
			const bool  &candidate_algorithm,
			const bool  &enable_outlier_rejection) : m_sigma_cutoff(sigma_cutoff),
								 m_OR_sigma_multiplier(OR_sigma_multiplier),
								 m_nHarmonics(nHarmonics),
								 m_export_powers(export_powers),
								 m_candidate_algorithm(candidate_algorithm),
								 m_enable_outlier_rejection(enable_outlier_rejection) {
      
    }
    
    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }

    float OR_sigma_multiplier() const {
      return m_OR_sigma_multiplier;
    }

    int nHarmonics() const {
      return m_nHarmonics;
    }

    int export_powers() const {
      return m_export_powers;
    }

    bool m_candidate_algorithm() const {
      return m_candidate_algorithm;
    }

    bool enable_outlier_rejection() const {
      return m_enable_outlier_rejection;
    }
    
  private:
    float m_sigma_cutoff;
    float m_OR_sigma_multiplier;
    int   m_nHarmonics;
    int   m_export_powers;
    bool  m_candidate_algorithm;
    bool  m_enable_outlier_rejection;
  };

}

#endif /* ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP */
