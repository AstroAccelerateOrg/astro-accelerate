//
//  aa_fdas_strategy.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 03/12/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP

#include <stdio.h>

#include "aa_strategy.hpp"
#include "aa_fdas_plan.hpp"

namespace astroaccelerate {

  /**
   * Class for aa_fdas_strategy, used to configure the fourier domain accelerated search (fdas).
   */

  class aa_fdas_strategy : public aa_strategy {
  public:
    aa_fdas_strategy() : m_sigma_cutoff(0.0), m_ready(false) {
      
    }
    
    aa_fdas_strategy(const aa_fdas_plan &fdas_plan) : m_sigma_cutoff(fdas_plan.sigma_cutoff()),
						      m_narrow(fdas_plan.narrow()),
						      m_wide(fdas_plan.wide()),
						      m_aggression(fdas_plan.aggression()),
						      m_num_boots(fdas_plan.num_boots()),
						      m_num_trial_bins(fdas_plan.num_trial_bins()),
						      m_navdms(fdas_plan.navdms()),
						      m_nsearch(fdas_plan.nsearch()),
						      m_ready(true) {
      
    }

    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }

    float narrow() const {
      return m_narrow;
    }

    float wide() const {
      return m_wide;
    }

    float aggression() const {
      return m_aggression;
    }

    int num_boots() const {
      return m_num_boots;
    }

    int num_trial_bins() const {
      return m_num_trial_bins;
    }

    int navdms() const {
      return m_navdms;
    }

    int	nsearch() const	{
      return m_nsearch;
    }

    bool setup() {
      return false;
    }

    bool ready() const {
      return m_ready;
    }
  private:
    float m_sigma_cutoff;
    float m_narrow;
    float m_wide;
    float m_aggression;
    int	  m_num_boots;
    int	  m_num_trial_bins;
    int	  m_navdms;
    int	  m_nsearch;

    bool m_ready;
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP
