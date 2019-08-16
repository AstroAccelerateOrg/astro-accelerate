#ifndef ASTRO_ACCELERATE_AA_ANALYSIS_PLAN_HPP
#define ASTRO_ACCELERATE_AA_ANALYSIS_PLAN_HPP

#include <stdio.h>

#include "aa_ddtr_strategy.hpp"

namespace astroaccelerate {
  /** \class aa_analysis_plan aa_analysis_plan.hpp "include/aa_analysis_plan.hpp"
   * \brief Class to set an aa_analysis_plan.
   * \details An analysis plan is required in order to create an analsyis strategy.   
   * \author AstroAccelerate.
   * \date 14 August 2019.
   */
  class aa_analysis_plan {
  public:

    enum class selectable_candidate_algorithm : int {
		peak_find = 0,
		threshold = 1,
		peak_filtering = 2
	};
    
    /** \brief Trivial constructor for aa_analysis_plan.
     * \details A trivial is used to return trivial plans which can only create strategies which cannot be ready.
     */
    aa_analysis_plan() {
      
    }
    
    /** \brief Constructor for aa_analysis_plan.
     * \details All parameters must be provided on construction and can only be set on construction.
     */
    aa_analysis_plan(const aa_ddtr_strategy &ddtr_strategy,
		     const float &sigma_cutoff,
		     const float &sigma_constant,
		     const float &max_boxcar_width_in_sec,
		     const selectable_candidate_algorithm &candidate_algorithm,
		     const bool &enable_msd_baseline_noise) : m_ddtr_strategy(ddtr_strategy),
							      m_sigma_cutoff(sigma_cutoff),
							      m_sigma_constant(sigma_constant),
							      m_max_boxcar_width_in_sec(max_boxcar_width_in_sec),
							      m_candidate_algorithm(candidate_algorithm),
							      m_enable_msd_baseline_noise(enable_msd_baseline_noise) {
      
    }
    
    /** \returns The aa_ddtr_strategy that was configured with an instance of aa_analysis_strategy. */
    const aa_ddtr_strategy ddtr_strategy() const {
      return m_ddtr_strategy;
    }

    /** \returns The sigma_cutoff that was configured with an instance of aa_analysis_strategy. */
    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }

    /** \returns The sigma_constant that was configured with an instance of aa_analysis_strategy. */
    float sigma_constant() const {
      return m_sigma_constant;
    }

    /** \returns The max_boxcar_width_in_sec that was configured with an instance of aa_analysis_strategy. */
    float max_boxcar_width_in_sec() const {
      return m_max_boxcar_width_in_sec;
    }

    /** \returns The setting for the selectable candidate algorithm (currently on or off) that was configured with an instance of aa_analysis_strategy. */
    aa_analysis_plan::selectable_candidate_algorithm candidate_algorithm() const {
      return m_candidate_algorithm;
    }

    /** \returns A boolean indicating whether the msd_baseline_noise algorithm will be enabled, for an instance of aa_analysis_strategy. */
    bool enable_msd_baseline_noise() const {
      return m_enable_msd_baseline_noise;
    }

  private:
    aa_ddtr_strategy    m_ddtr_strategy; /**< The instance of aa_ddtr_strategy configured for an instance off aa_analysis_plan. */
    float               m_sigma_cutoff; /**< The user defined sigma_cutoff. */
    float               m_sigma_constant; /**< The user defined sigma_constnat. */
    float               m_max_boxcar_width_in_sec; /**< The user defined max_boxcar_width_in_sec. */
    aa_analysis_plan::selectable_candidate_algorithm m_candidate_algorithm; /**< The user defined setting for the selectable_candidate_algorithm (currently on or off). */
    bool                m_enable_msd_baseline_noise; /**< The user defined boolean setting to enable/disable the msd_baseline_noise reduction algorithm. */
  };
  
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_ANALYSIS_PLAN_HPP
