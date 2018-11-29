//
//  aa_periodicity_strategy.hpp
//  aapipeline
//
//  Created by Cees Carels on Tuesday 23/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_PERIODICITY_STRATEGY_HPP
#define ASTRO_ACCELERATE_PERIODICITY_STRATEGY_HPP

#include <stdio.h>

#include "aa_strategy.hpp"
#include "aa_periodicity_plan.hpp"

namespace astroaccelerate {
  
  class aa_periodicity_strategy : public aa_strategy {
  public:
    aa_periodicity_strategy() : m_sigma_cutoff(0),
				m_OR_sigma_multiplier(0),
				m_nHarmonics(0),
				m_export_powers(0),
				m_candidate_algorithm(false),
				m_enable_outlier_rejection(false),
				m_ready(false) {
      
    }
    
    aa_periodicity_strategy(const aa_periodicity_plan &plan) : m_sigma_cutoff(plan.sigma_cutoff()),
							       m_OR_sigma_multiplier(plan.OR_sigma_multiplier()),
							       m_nHarmonics(plan.nHarmonics()),
							       m_export_powers(plan.export_powers()),
							       m_candidate_algorithm(plan.candidate_algorithm()),
							       m_enable_outlier_rejection(plan.enable_outlier_rejection()),
							       m_ready(false) {
      if((m_nHarmonics > 0) && (m_OR_sigma_multiplier > 0) && (m_export_powers > 0)) {
	m_ready = true;
      }      
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

    bool candidate_algorithm() const {
      return m_candidate_algorithm;
    }
  
    bool enable_outlier_rejection() const {
      return m_enable_outlier_rejection;
    }

    bool ready() const {
      return m_ready;
    }

    bool setup() {
      return m_ready;
    }

    bool print_parameters() const {
      printf("Periodicity - sigma_cutoff %f\n", m_sigma_cutoff);
      printf("Periodicity - OR_sigma_multiplier %f\n", m_OR_sigma_multiplier);
      printf("Periodicity - nHarmonics %d\n", m_nHarmonics);
      printf("Periodicity - export_powers %d\n", m_export_powers);
      printf("Periodicity - candidate_algorithm %d\n", m_candidate_algorithm);
      printf("Periodicity - enable_outlier_rejection %d\n", m_enable_outlier_rejection);
      return true;
    }
    
  private:
    float m_sigma_cutoff;
    float m_OR_sigma_multiplier;
    int   m_nHarmonics;
    int   m_export_powers;
    bool  m_candidate_algorithm;
    bool  m_enable_outlier_rejection;

    bool m_ready;
  };

} //namespace astroaccelerate

#endif /* ASTRO_ACCELERATE_PERIODICITY_STRATEGY_HPP */
