#ifndef ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP

#include <stdio.h>

#include "aa_strategy.hpp"
#include "aa_fdas_plan.hpp"

namespace astroaccelerate {

  /**
   * \class aa_fdas_strategy aa_fdas_strategy.hpp "include/aa_fdas_strategy.hpp" 
   * \brief Class that receives an aa_fdas_plan object, and produces an aa_fdas_strategy object.
   * \details The strategy is calculated based on the plan.
   * \details An fdas strategy is required for any pipeline running the fdas component.
   * \author Cees Carels.
   * \date 3 December 2018.
   */

  class aa_fdas_strategy : public aa_strategy {
  public:
    /** Trivial constructor for aa_fdas_strategy. */
    aa_fdas_strategy() : m_sigma_cutoff(0.0),
			 m_sigma_constant(0.0),
			 m_narrow(0.0),
			 m_wide(0.0),
			 m_aggression(0.0),
			 m_num_boots(0),
			 m_num_trial_bins(0),
			 m_navdms(0),
			 m_enable_msd_baseline_noise(false),
			 m_ready(false) {
      
    }
    
    /** Constructor aa_fdas_strategy that initialises al member variables. */
    aa_fdas_strategy(const aa_fdas_plan &fdas_plan) : m_sigma_cutoff(fdas_plan.sigma_cutoff()),
						      m_sigma_constant(fdas_plan.sigma_constant()),
						      m_narrow(fdas_plan.narrow()),
						      m_wide(fdas_plan.wide()),
						      m_aggression(fdas_plan.aggression()),
						      m_num_boots(fdas_plan.num_boots()),
						      m_num_trial_bins(fdas_plan.num_trial_bins()),
						      m_navdms(fdas_plan.navdms()),
						      m_enable_msd_baseline_noise(fdas_plan.enable_msd_baseline_noise()),
						      m_ready(false) {

      /** Parse input. Invalid input means the ready state will not be set to true. */
      if((m_sigma_cutoff > 0)
	 && (m_sigma_constant > 0)
	 && (m_narrow > 0)
	 && (m_wide > 0)
	 && (m_aggression > 0)
	 && (m_num_boots > 0)
	 && (m_num_trial_bins > 0)
	 && (m_navdms > 0)) {
	m_ready = true;
      }
    }

    /** \returns The name of this mdoule. */
    std::string name() const {
      return "fdas_strategy";
    }
    
    /** \returns The strategy determined sigma_cutoff. */
    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }

    /** \returns The strategy determined sigma_constant. */
    float sigma_constant() const {
      return m_sigma_constant;
    }

    /** \returns The strategy determined narrow setting. */
    float narrow() const {
      return m_narrow;
    }

    /** \returns The strategy determined wide setting. */
    float wide() const {
      return m_wide;
    }

    /** \returns The strategy determined aggression setting. */
    float aggression() const {
      return m_aggression;
    }

    /** \returns The strategy determined num_boots setting. */
    int num_boots() const {
      return m_num_boots;
    }

    /** \returns The strategy determined num_trial_bins setting. */
    int num_trial_bins() const {
      return m_num_trial_bins;
    }

    /** \returns The strategy determined navdms setting. */
    int navdms() const {
      return m_navdms;
    }
    
    /**
     * \returns an integer to indicate whether the msd baseline noise reduction algorithm will be enabled or disabled. 0 for off (false), 1 for on (true).
     * \details At the moment, this setting has no effect.
     **/
    int enable_msd_baseline_noise() const {
      return (m_enable_msd_baseline_noise) ? 1 : 0;
    }
    
    /** \brief Performs any setup still needed for the strategy.
     * \returns A boolean indicating whether the setup was successful.
     */
    bool setup() {
      return ready();
    }

    /** \returns The ready state of the instance of the fdas_strategy. */
    bool ready() const {
      return m_ready;
    }

    /**
     * \brief Static member function that prints member variables for a provided aa_fdas_strategy instance.
     * \returns A boolean to indicate whether the printing was successful.
     */
    static bool print_info(const aa_fdas_strategy &strategy) {
      LOG(log_level::dev_debug, "FDAS STRATEGY INFORMATION");
      LOG(log_level::dev_debug, "fdas sigma_cutoff:\t\t\t" + std::to_string(strategy.sigma_cutoff()));
      LOG(log_level::dev_debug, "fdas_sigma_constant:\t\t\t" + std::to_string(strategy.sigma_constant()));
      LOG(log_level::dev_debug, "fdas narrow:\t\t\t\t" + std::to_string(strategy.narrow()));
      LOG(log_level::dev_debug, "fdas wide:\t\t\t\t" + std::to_string(strategy.wide()));
      LOG(log_level::dev_debug, "fdas aggression:\t\t\t\t" + std::to_string(strategy.aggression()));
      LOG(log_level::dev_debug, "fdas num_boots:\t\t\t\t" + std::to_string(strategy.num_boots()));
      LOG(log_level::dev_debug, "fdas num_trial_bins:\t\t\t" + std::to_string(strategy.num_trial_bins()));
      LOG(log_level::dev_debug, "fdas navdms:\t\t\t\t" + std::to_string(strategy.navdms()));
      return true;
    }
  private:
    float m_sigma_cutoff; /**< The strategy determined sigma_cutoff setting. */
    float m_sigma_constant; /**< The strategy determined sigma_constant setting. */
    float m_narrow; /**< The strategy determined narrow setting. */
    float m_wide; /**< The strategy determined wide setting. */
    float m_aggression; /**< The strategy determined aggression setting. */
    int	  m_num_boots; /** The strategy determined num_boots setting. */
    int	  m_num_trial_bins; /** The strategy determined num_trial_bins setting. */
    int	  m_navdms; /** The strategy determined navdms setting. */
    bool  m_enable_msd_baseline_noise; /** Flag for enabling/disabling msd_baseline_noise reduction algorithm. */
    
    bool m_ready; /** The ready state of the instance of the fdas_strategy. */
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP
