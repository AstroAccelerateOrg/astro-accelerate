#ifndef ASTRO_ACCELERATE_AA_DDTR_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_DDTR_STRATEGY_HPP

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "aa_strategy.hpp"

#include "aa_params.hpp"
#include "aa_ddtr_plan.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_device_info.hpp"

#include "aa_log.hpp"

namespace astroaccelerate {

  /**
   * \class aa_ddtr_strategy aa_ddtr_strategy.hpp "include/aa_ddtr_strategy.hpp"
   * \brief Class that receives an aa_ddtr_plan object, and produces an aa_ddtr_strategy object.
   * \details The strategy is calculated based on the plan.
   * \details A ddtr strategy is required for any pipeline running the dedispersion component.
   * \details The code that does this was formerly known as "stratagy" in "host_stratagy.cu".
   * \details Further implementation details are provided in aa_ddtr_strategy.cpp.
   * \warning Upon construction, the user must specify whether aa_pipeline::component::analysis will be run in the pipeline.
   * \warning Incorrectly setting enable_analysis upon construction may lead to inefficient memory usage, or undefined behaviour.
   * \author AstroAccelerate team.
   * \date 23 October 2018.
   */
  class aa_ddtr_strategy : public aa_strategy {
  public:
    /** \brief Implementation of trivial constructor in source file. */
    aa_ddtr_strategy();

    /** \brief Implementation of non-trvial constructor in source file. All parameters are set once on construction. */
    aa_ddtr_strategy(const aa_ddtr_plan &plan, const aa_filterbank_metadata &metadata, const size_t &free_memory, const bool &enable_analysis, aa_device_info *selected_device);

    /** \brief Destructor for aa_ddtr_strategy. */
    virtual ~aa_ddtr_strategy() {

    }

    /** \returns The name of the module. */
    std::string name() const {
      return "ddtr_strategy";
    }
    
    /** Implementation in source file. */
    bool setup();
    
    /** \returns maxshift as set by strategy. */
    int maxshift() const {
      return m_maxshift;
    }
    
    /** \returns t_processed as set by strategy. */
    const std::vector<std::vector<int>>& t_processed() const {
      return m_t_processed;
    }
	
	size_t nProcessedTimesamples() {
		size_t sum = 0;
		if(m_t_processed.size()>0){
			for(int f=0; f<(int)m_t_processed[0].size();f++){
				sum = sum + m_t_processed[0][f];
			}
		}
		return(sum);
	}
		
    
    /** \returns dmshifts as set by strategy. */
    std::vector<float> dmshifts() const {
      return m_dmshifts;
    }
	
    /** \returns size of the custom bandpass normalization. */
    size_t bandpass_normalization_size() const {
      return m_bandpass_normalization.size();
    }
	
    /** \returns size of the custom bandpass normalization. */
    const float* bandpass_normalization_pointer() const {
      return m_bandpass_normalization.data();
    }
	
    /** \returns custom bandpass normalization as set by the user in ddtr plan. */
    std::vector<float> bandpass_normalization() const {
      return m_bandpass_normalization;
    }
    
    /** \returns dm at an index in a std::vector of dm. */
    const aa_ddtr_plan::dm dm(const size_t &i) const {
      return str_dm.at(i);
    }
    
    /** \returns the number of dm ranges. */
    size_t get_nRanges() const {
      return str_dm.size();
    }

    /** \returns The number of dedispersion measures (dm). */
    size_t ndms_size() const {
      return m_ndms.size();
    }
    
    /** \returns The number of dm at an index in a std::vector of dm. */
    int ndms(const size_t &i) const {
      return m_ndms.at(i);
    }
    
    /** \returns The maximum dm. */
    int max_ndms() const {
      return m_max_ndms;
    }
    
    /** \returns A pointer to the ndms data. */
    const int* ndms_data() const {
      return m_ndms.data();
    }

    /** return the DM low for specified range */
    int dm_low(const int range) const{
	    return str_dm.at(range).low;
    }
    
    /** \returns The number of time chunks. */
    int num_tchunks() const {
      return m_num_tchunks;
    }

    /** \returns The strategy configured power. */
    float power() const {
      return m_power;
    }

    /**
     * \returns an integer to indicate whether the msd baseline noise reduction algorithm will be enabled or disabled. 0 for off (false), 1 for on (true).
     * \details At the moment, this setting has no effect.
     */
    int enable_msd_baseline_noise() const {
      return (m_enable_msd_baseline_noise) ? 1 : 0;
    }
    
    /** \returns The boolean ready state of the strategy instance (true for ready, false for not ready). */
    bool ready() const {
      return m_ready;
    }

    /** \returns The aa_filterbank_metadata instance bound to the instance of strategy. */
    const aa_filterbank_metadata metadata() const {
      return m_metadata;
    }

    /** \returns A boolean indicating whether the ddtr strategy instance was configured to run aa_pipeline::component::analysis. */
    bool configured_for_analysis() const {
      return m_configured_for_analysis;
    }

    /** \brief Static member function that prints member data for an aa_ddtr_strategy object. */
    static bool print_info(const aa_ddtr_strategy &strategy) {
      LOG(log_level::dev_debug, "DDTR STRATEGY INFORMATION:");
      LOG(log_level::dev_debug, "Optimization settings:");
      LOG(log_level::dev_debug, "\t\tUNROLL:\t\t\t" + std::to_string(UNROLLS));
      LOG(log_level::dev_debug, "\t\tSNUMREG:\t\t" + std::to_string(SNUMREG));
      LOG(log_level::dev_debug, "\t\tSDIVINT:\t\t" + std::to_string(SDIVINT));
      LOG(log_level::dev_debug, "\t\tSDIVINDM:\t\t" + std::to_string(SDIVINDM));
      LOG(log_level::dev_debug, "ddtr+analysis:\t\t" +  (strategy.configured_for_analysis() ? std::string("true") : std::string("false")));
      LOG(log_level::dev_debug, "ddtr dm ranges:\t\t" + std::to_string(strategy.get_nRanges()));
      for(size_t i = 0; i < strategy.get_nRanges(); i++) {
        const aa_ddtr_plan::dm tmp = strategy.dm(i);
        LOG(log_level::dev_debug, "     dm (low,high,step,inBin,outBin) " +
	      std::to_string(tmp.low) + "," + std::to_string(tmp.high) + "," + std::to_string(tmp.step)
	      + "," + std::to_string(tmp.inBin) + "," + std::to_string(tmp.outBin)
	    );
      }

      LOG(log_level::dev_debug, "ddtr max_ndms:\t\t" + std::to_string(strategy.max_ndms()));
      LOG(log_level::dev_debug, "ddtr ndms elements:");
      for(size_t i = 0; i < strategy.ndms_size(); i++) {
        LOG(log_level::dev_debug, "     ndms[" + std::to_string(i) + "]:\t\t" + std::to_string(strategy.ndms(i)));
      }
      
      LOG(log_level::dev_debug, "ddtr maxshift:\t\t" + std::to_string(strategy.maxshift()));
      LOG(log_level::dev_debug, "ddtr num_tchunks:\t" + std::to_string(strategy.num_tchunks()));
      LOG(log_level::dev_debug, "ddtr max_ndms:\t\t" + std::to_string(strategy.max_ndms()));
      LOG(log_level::dev_debug, "t_processed size:\t" + std::to_string(strategy.t_processed().size()));
      LOG(log_level::dev_debug, "t_processed elements:");
      for(size_t i = 0; i < strategy.t_processed().size(); i++) {
        for(size_t j = 0; j < strategy.t_processed().at(i).size(); j++) {
          LOG(log_level::dev_debug, "     t_processed[" + std::to_string(i) +"][" + std::to_string(j) + "]:\t" + std::to_string(strategy.t_processed()[i][j]));
        }
      }
      LOG(log_level::dev_debug, "power:\t\t\t" + std::to_string(strategy.power()));

      return true;
    }
    
  private:
    bool strategy(const aa_ddtr_plan &plan, const size_t &free_memory, const bool &enable_analysis);
    bool m_ready; /**< The ready state of the ddtr strategy. */
    bool m_strategy_already_calculated; /**< A flag to indicate whether the strategy for the instance has already been allocated. */
	aa_device_info *m_selected_device;

    bool m_configured_for_analysis; /**< A flag to indicate whether the strategy will be configured to run aa_pipeline::component::analysis. */
    bool is_setup;  /**< A flag to indicate whether the setup method has been called. */

    aa_filterbank_metadata m_metadata;
    std::vector<int> m_ndms;
    std::vector<float> m_dmshifts;
    std::vector<aa_ddtr_plan::dm> str_dm;
    int m_maxshift;       /**< Is used for assignment and assigning. */
    int m_num_tchunks;    /**< Is used for assignment only. */
    int m_total_ndms;     /**< Is used for assignment only. */
    float m_max_dm;       /**< Is used for assignment only. */
    int m_maxshift_high;  /**< Is used for assignment and assigning in this method. */
    
    int m_max_ndms;       /**< This variable is set to 0 in main.cpp and never used until here. */
    float m_power;        /**< Strategy determined power from the user set aa_ddtr_plan instance of power. */
    std::vector<std::vector<int>> m_t_processed; /**< Is allocated in this class, and used elsewhere in the pipeline. */
    float ***output_buffer; /**< \brief 3D array that contains the output. \deprecated Has been moved to permitted_pipeline classes. Remove from source and header files.*/
    bool m_enable_msd_baseline_noise; /** Flag that enables or disables the use of msd baseline noise. */
	
	std::vector<float> m_bandpass_normalization;
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DDTR_STRATEGY_HPP
