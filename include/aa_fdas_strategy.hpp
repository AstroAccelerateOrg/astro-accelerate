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
   * \author AstroAccelerate.
   * \date 1 July 2019.
   */

  class aa_fdas_strategy : public aa_strategy {
  private:
    float m_sigma_cutoff; /**< The strategy determined sigma_cutoff setting. */
    float m_sigma_constant; /**< The strategy determined sigma_constant setting. */
    bool  m_enable_msd_baseline_noise; /** Flag for enabling/disabling msd_baseline_noise reduction algorithm. */
    
    bool m_ready; /** The ready state of the instance of the fdas_strategy. */
  public:
    /** Trivial constructor for aa_fdas_strategy. */
    aa_fdas_strategy() : m_sigma_cutoff(0.0),
			 m_sigma_constant(0.0),
			 m_enable_msd_baseline_noise(false),
			 m_ready(false) {
      
    }
    
    /** Constructor aa_fdas_strategy that initialises al member variables. */
    aa_fdas_strategy(const aa_fdas_plan &fdas_plan) : m_sigma_cutoff(fdas_plan.sigma_cutoff()),
						      m_sigma_constant(fdas_plan.sigma_constant()),
						      m_enable_msd_baseline_noise(fdas_plan.enable_msd_baseline_noise()),
						      m_ready(false) {

      /** Parse input. Invalid input means the ready state will not be set to true. */
		if((m_sigma_cutoff > 0) && (m_sigma_constant > 0)) {
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
      LOG(log_level::dev_debug, "\t\tNEXP:\t\t\t" + std::to_string(NEXP));
      LOG(log_level::dev_debug, "\t\tPOTWO:\t\t\t" + std::to_string(POTWO));
      LOG(log_level::dev_debug, "\t\tKERNLEN:\t\t" + std::to_string(KERNLEN));
      LOG(log_level::dev_debug, "\t\tACCEL_STEP:\t\t" + std::to_string(ACCEL_STEP));
      LOG(log_level::dev_debug, "\t\tACCEL_STEP_R:\t\t" + std::to_string(ACCEL_STEP_R));
      LOG(log_level::dev_debug, "\t\tZMAX:\t\t\t" + std::to_string(ZMAX));
      LOG(log_level::dev_debug, "\t\tNKERN:\t\t\t" + std::to_string(NKERN));
      return true;
    }
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP
