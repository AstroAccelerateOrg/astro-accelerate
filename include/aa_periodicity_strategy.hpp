#ifndef ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP

#include <stdio.h>

#include "aa_strategy.hpp"
#include "aa_periodicity_plan.hpp"
#include "aa_log.hpp"

namespace astroaccelerate {

  /**
   * \class aa_periodicity_strategy aa_periodicity_strategy.hpp "include/aa_periodicity_strategy.hpp"
   * \brief Class that receives an aa_periodicity_plan object, and produces an aa_periodicity_strategy object.
   * \details A periodicity strategy is required for any pipeline running the periodicity component.
   * \author Cees Carels.
   * \date 23 October 2018.
   */
  
  class aa_periodicity_strategy : public aa_strategy {
  public:

    /** \brief Trivial constructor for aa_periodicity_strategy, which can never have a ready state equal to true. */
    aa_periodicity_strategy() : m_sigma_cutoff(0),
				m_sigma_constant(0),
				m_nHarmonics(0),
				m_export_powers(0),
				m_candidate_algorithm(false),
				m_enable_msd_baseline_noise(false),
				m_ready(false) {
      
    }
    
    /** \brief Constructor for aa_periodicity_strategy that sets all member variables upon construction. */
    aa_periodicity_strategy(const aa_periodicity_plan &plan) : m_sigma_cutoff(plan.sigma_cutoff()),
							       m_sigma_constant(plan.sigma_constant()),
							       m_nHarmonics(plan.nHarmonics()),
							       m_export_powers(plan.export_powers()),
							       m_candidate_algorithm(plan.candidate_algorithm()),
							       m_enable_msd_baseline_noise(plan.enable_msd_baseline_noise()),
							       m_ready(false) {
      /** Parse user input, if the user input is not valid, then the ready state will not become true. */
      if((m_nHarmonics > 0) && (m_sigma_constant > 0) && (m_export_powers >= 0)) {
	m_ready = true;
      }
      else {
	LOG(log_level::warning, "Invalid periodicity strategy parameters. Check the aa_periodicity_plan input parameters.");
	print_info(*this);
      }
    }

    /** \returns The name of the module. */
    std::string name() const {
      return "periodicity_strategy";
    }
    
    /** \returns The strategy determined sigma_cutoff. */
    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }

    /** \returns The strategy determined sigma_constant. */
    float sigma_constant() const {
      return m_sigma_constant;
    }

    /** \returns The strategy determined nHarmonics. */
    int nHarmonics() const {
      return m_nHarmonics;
    }

    /** \returns The strategy determined export_powers. */
    int export_powers() const {
      return m_export_powers;
    }

    /** \returns The strategy determined candidate algorithm flag. */
    bool candidate_algorithm() const {
      return m_candidate_algorithm;
    }
  
    /** \returns an integer to indicate whether the msd baseline noise reduction algorithm will be enabled or disabled. 0 for off (false), 1 for on (true). */
    int enable_msd_baseline_noise() const {
      return (m_enable_msd_baseline_noise) ? 1 : 0;
    }
    
    /** \returns The ready state of the instance of the class. */
    bool ready() const {
      return m_ready;
    }
    
    /** \brief Performs any remaining setup needed. 
     * \returns Whether the operation was successful or not, at the moment this is the ready state of the instance. */
    bool setup() {
      return ready();
    }

    /** \brief Print the member data of the instnace.
     * \returns A boolean flag to indicate whether the operation completed successfully (true) or unsuccessfully (false). */
    bool print_parameters() const {
      printf("Periodicity - sigma_cutoff %f\n", m_sigma_cutoff);
      printf("Periodicity - sigma_constant %f\n", m_sigma_constant);
      printf("Periodicity - nHarmonics %d\n", m_nHarmonics);
      printf("Periodicity - export_powers %d\n", m_export_powers);
      printf("Periodicity - candidate_algorithm %d\n", m_candidate_algorithm);
      printf("Periodicity - enable_msd_baseline_noise %d\n", m_enable_msd_baseline_noise);
      return true;
    }

    /** Static member function that prints member variables for a provided aa_periodicity_strategy. */
    static bool print_info(const aa_periodicity_strategy &strategy) {
      LOG(log_level::dev_debug, "PERIODICITY STRATEGY INFORMATION:");
      LOG(log_level::dev_debug, "periodicity sigma_cutoff:\t\t" + std::to_string(strategy.sigma_cutoff()));
      LOG(log_level::dev_debug, "periodicity sigma_constant:\t\t" + std::to_string(strategy.sigma_constant()));
      LOG(log_level::dev_debug, "periodicity nHarmonics:\t\t\t" + std::to_string(strategy.nHarmonics()));
      LOG(log_level::dev_debug, "periodicity export_powers:\t\t" + std::to_string(strategy.export_powers()));
      LOG(log_level::dev_debug, "periodicity candidate_algorithm:\t\t" + (strategy.candidate_algorithm() ? std::string("true") : std::string("false")));
      LOG(log_level::dev_debug, "periodicity enable_msd_baseline_noise:\t" + (strategy.enable_msd_baseline_noise() ? std::string("true") : std::string("false")));
      return true;
    }
    
  private:
    float m_sigma_cutoff; /** Strategy determined sigma_cutoff. */
    float m_sigma_constant; /** Strategy determined sigma_constant. */
    int   m_nHarmonics; /** Strategy determined nHarmonics. */
    int   m_export_powers; /** Strategy determined export_powers. */
    bool  m_candidate_algorithm; /** Strategy determined candidate algorithm flag. */
    bool  m_enable_msd_baseline_noise; /** Strategy determined msd baseline noise flag. */

    bool m_ready; /** Ready state of the instance. */
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP
