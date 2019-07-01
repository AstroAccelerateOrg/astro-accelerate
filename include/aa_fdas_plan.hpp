#ifndef ASTRO_ACCELERATE_AA_FDAS_PLAN_HPP
#define ASTRO_ACCELERATE_AA_FDAS_PLAN_HPP

namespace astroaccelerate {

  /**
   * \class aa_fdas_plan aa_fdas_plan.hpp "include/aa_fdas_plan.hpp"  
   * \brief Class for aa_fdas_plan, used to configure the fourier domain accelerated search (fdas).
   * \details An fdas plan is required to create an fdas strategy.
   * \author AstroAccelerate.
   * \date 3 December 2018.
   */

  class aa_fdas_plan {
  public:
    /**
     * Trivial constructor for aa_fdas_plan.
     */
    aa_fdas_plan() {
      
    }

    /**
     * Constructor for aa_fdas_plan that initialises all member variables.
     */
    aa_fdas_plan(const float &sigma_cutoff,
		 const float &sigma_constant,
		 const bool  &enable_msd_baseline_noise) : m_sigma_cutoff(sigma_cutoff),
							   m_sigma_constant(sigma_constant),
							   m_enable_msd_baseline_noise(enable_msd_baseline_noise) {
      
    }

    /** \returns User selected sigma_cutoff. */
    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }

    /** \returns User selected sigma_constant. */
    float sigma_constant() const {
      return m_sigma_constant;
    }

    /**
     * \returns A boolean indicating whether the msd_baseline_noise algorithm will be enabled, for an instance of aa_analysis_strategy.
     * \details At the moment, this setting has no effect.
     */
    bool enable_msd_baseline_noise() const {
      return m_enable_msd_baseline_noise;
    }
    
  private:
    float m_sigma_cutoff; /**< User selected sigma_cutoff. */
    float m_sigma_constant; /**< User selected sigma_constant. */
    bool  m_enable_msd_baseline_noise; /**< User selected flag for enabling / disabling msd_baseline_noise. */
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_FDAS_PLAN_HPP
