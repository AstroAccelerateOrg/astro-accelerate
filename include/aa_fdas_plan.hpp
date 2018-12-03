//
//  aa_fdas_plan.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 03/12/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_AA_FDAS_PLAN_HPP
#define ASTRO_ACCELERATE_AA_FDAS_PLAN_HPP

namespace astroaccelerate {

  /**
   * Class for aa_fdas_plan, used to configure the fourier domain accelerated search (fdas).
   */

  class aa_fdas_plan {
  public:
    aa_fdas_plan() {
      
    }

    aa_fdas_plan(const float &sigma_cutoff,
		 const int &num_boots,
		 const int &num_trial_bins,
		 const int &navdms,
		 const float &narrow,
		 const float &wide,
		 const int &nsearch,
		 const float &aggression) : m_sigma_cutoff(sigma_cutoff),
					    m_narrow(narrow),
					    m_wide(wide),
					    m_aggression(aggression),
					    m_num_boots(num_boots),
					    m_num_trial_bins(num_trial_bins),
					    m_navdms(navdms),
					    m_nsearch(nsearch) {
      
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

    int nsearch() const {
      return m_nsearch;
    }
  private:
    float m_sigma_cutoff;
    float m_narrow;
    float m_wide;
    float m_aggression;
    int   m_num_boots;
    int   m_num_trial_bins;
    int   m_navdms;
    int   m_nsearch;
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_FDAS_PLAN_HPP
