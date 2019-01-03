#ifndef ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP
#define ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP

#include <stdio.h>

namespace astroaccelerate {

  /**
   * \class aa_periodicity_plan aa_periodicity_plan.hpp "include/aa_periodicity_plan.hpp"
   * \brief Class to set a periodicity plan.
   * \details A periodicity plan is required to create a periodicity strategy.
   * \author Cees Carels.
   * \date 23 October 2018.
   */

  class aa_periodicity_plan {
  public:

    /** \brief Trivial constructor for aa_periodicity_plan. */
    aa_periodicity_plan() : m_sigma_cutoff(0),
			    m_OR_sigma_multiplier(0),
			    m_nHarmonics(0),
			    m_export_powers(0),
			    m_candidate_algorithm(false),
			    m_enable_outlier_rejection(false) {
      
    }
    
    /** \brief Constructor for aa_periodicity_plan that sets all member data on construction. */
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
    
    /** \returns User selected sigma_cutoff. */
    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }

    /** \returns User selected OR_sigma_multiplier. */
    float OR_sigma_multiplier() const {
      return m_OR_sigma_multiplier;
    }

    /** \returns User selected nHarmonics. */
    int nHarmonics() const {
      return m_nHarmonics;
    }

    /** \returns User selected export_powers. */
    int export_powers() const {
      return m_export_powers;
    }

    /** \returns A boolean flag to indicate whether the candidate algorithm was selected. True indicates that it is, false indicates that it is not. */
    bool candidate_algorithm() const {
      return m_candidate_algorithm;
    }
    
    /** \returns A boolean flag to indicate whether the outlier rejection is selected. True indicates that it is, false indicates that it is not. */
    bool enable_outlier_rejection() const {
      return m_enable_outlier_rejection;
    }
    
  private:
    float m_sigma_cutoff; /** User selected sigma_cutoff. */
    float m_OR_sigma_multiplier; /** User selected OR_sigma_multiplier. */
    int   m_nHarmonics; /** User selected nHarmonics. */
    int   m_export_powers; /** User selected export_powers. */
    bool  m_candidate_algorithm; /** User selected flag to enable or disable the use of the candidate_algorithm. */
    bool  m_enable_outlier_rejection; /** User selected flag to enable or disable outlier rejection. */
  };
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP
