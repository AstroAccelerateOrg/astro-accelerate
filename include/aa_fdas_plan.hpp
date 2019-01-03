#ifndef ASTRO_ACCELERATE_AA_FDAS_PLAN_HPP
#define ASTRO_ACCELERATE_AA_FDAS_PLAN_HPP

namespace astroaccelerate {

  /**
   * \class aa_fdas_plan aa_fdas_plan.hpp "include/aa_fdas_plan.hpp"  
   * \brief Class for aa_fdas_plan, used to configure the fourier domain accelerated search (fdas).
   * \details An fdas plan is required to create an fdas strategy.
   * \author Cees Carels.
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

    /** \returns User selected sigma_cutoff. */
    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }

    /** \returns User selected narrow setting. */
    float narrow() const {
      return m_narrow;
    }

    /** \returns User selected wide setting. */
    float wide() const {
      return m_wide;
    }

    /** \returns User selected aggression setting. */
    float aggression() const {
      return m_aggression;
    }

    /** \returns User selected num_boots. */
    int num_boots() const {
      return m_num_boots;
    }

    /** \returns User selected number of trial bins. */
    int num_trial_bins() const {
      return m_num_trial_bins;
    }

    /** \returns User selected number of dms. */
    int navdms() const {
      return m_navdms;
    }

    /** \returns User selected nsearch parameter. */
    int nsearch() const {
      return m_nsearch;
    }
  private:
    float m_sigma_cutoff; /**< User selected sigma_cutoff. */
    float m_narrow; /**< User selected narrow setting. */
    float m_wide; /**< User selected wide setting. */
    float m_aggression; /**< User selected aggression setting. */
    int   m_num_boots; /**< User selected num_boots. */
    int   m_num_trial_bins; /**< User selected num_trial_bins. */
    int   m_navdms; /**< User selected dms. */
    int   m_nsearch; /**< User selected nsearch parameter. */
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_FDAS_PLAN_HPP
