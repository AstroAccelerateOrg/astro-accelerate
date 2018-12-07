//
//  aa_ddtr_strategy.hpp
//  aapipeline
//
//  Created by Cees Carels on Tuesday 23/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_DDTR_STRATEGY_HPP
#define ASTRO_ACCELERATE_DDTR_STRATEGY_HPP

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "aa_strategy.hpp"

#include "params.hpp"
#include "aa_ddtr_plan.hpp"
#include "aa_filterbank_metadata.hpp"

/**
 * Class that receives an aa_ddtr_plan object, and produces an aa_ddtr_strategy object.
 * The strategy is calculated based on the plan.
 * The code that does this was formerly known as "stratagy" in "host_stratagy.cu".
 */

namespace astroaccelerate {

  class aa_ddtr_strategy : public aa_strategy {
  public:
    aa_ddtr_strategy();
    aa_ddtr_strategy(const aa_ddtr_plan &plan, const aa_filterbank_metadata &metadata, const size_t &free_memory, const bool &enable_analysis);
    ~aa_ddtr_strategy() {

    }

    std::string name() const {
      return "ddtr_strategy";
    }
    
    bool setup();
    
    int maxshift() const {
      return m_maxshift;
    }
    
    const std::vector<std::vector<int>>& t_processed() const {
      return m_t_processed;
    }
    
    std::vector<float> dmshifts() const {
      return m_dmshifts;
    }
    
    const aa_ddtr_plan::dm dm(const size_t &i) const {
      return str_dm.at(i);
    }
    
    size_t range() const {
      return str_dm.size();
    }

    size_t ndms_size() const {
      return m_ndms.size();
    }
    
    int ndms(const size_t &i) const {
      return m_ndms.at(i);
    }
    
    int max_ndms() const {
      return m_max_ndms;
    }
    
    const int* ndms_data() const {
      return m_ndms.data();
    }
    
    int num_tchunks() const {
      return m_num_tchunks;
    }
    
    bool ready() const {
      return m_ready;
    }

    const aa_filterbank_metadata metadata() const {
      return m_metadata;
    }

    bool configured_for_analysis() const {
      return m_configured_for_analysis;
    }

    static bool print_info(const aa_ddtr_strategy &strategy) {
      std::cout << "DDTR STRATEGY INFORMATION:" << std::endl;
      std::cout << "ddtr+analysis:\t\t" << (strategy.configured_for_analysis() ? "true" : "false") << std::endl;
      std::cout << "ddtr dm ranges:\t\t" << strategy.range() << std::endl;
      for(size_t i = 0; i < strategy.range(); i++) {
	const aa_ddtr_plan::dm tmp = strategy.dm(i);
	std::cout << "     dm (low,high,step,inBin,outBin) "
		  << tmp.low << "," << tmp.high << "," << tmp.step << ","
		  << tmp.inBin << "," << tmp.outBin << std::endl;
      }

      std::cout << "ddtr max_ndms:\t\t" << strategy.max_ndms() << std::endl;
      std::cout << "ddtr ndms elements:" << std::endl;
      for(size_t i = 0; i < strategy.ndms_size(); i++) {
	std::cout << "     ndms[" << i << "]:\t\t" << strategy.ndms(i) << std::endl;
      }
      
      std::cout << "ddtr maxshift:\t\t" << strategy.maxshift() << std::endl;
      std::cout << "ddtr num_tchunks:\t" << strategy.num_tchunks() << std::endl;
      std::cout << "ddtr max_ndms:\t\t" << strategy.max_ndms() << std::endl;
      std::cout << "t_processed size:\t" << strategy.t_processed().size() << std::endl;
      std::cout << "t_processed elements:" << std::endl;
      for(size_t i = 0; i < strategy.t_processed().size(); i++) {
	for(size_t j = 0; j < strategy.t_processed().at(i).size(); j++) {
	  std::cout << "     t_processed[" << i << "][" << j << "]:\t" << strategy.t_processed()[i][j] << std::endl;
	}
      }

      return true;
    }
    
  private:
    bool strategy(const aa_ddtr_plan &plan, const size_t &free_memory, const bool &enable_analysis);
    bool m_ready;
    bool m_strategy_already_calculated;

    bool m_configured_for_analysis;
    bool is_setup;  //Has setup been called already?

    aa_filterbank_metadata m_metadata;
    std::vector<int> m_ndms;
    std::vector<float> m_dmshifts;
    std::vector<aa_ddtr_plan::dm> str_dm;
    int m_maxshift;       //Is used for assignment and assigning
    int m_num_tchunks;    //Is used for assignment only
    int m_total_ndms;     //Is used for assignment only
    float m_max_dm;       //Is used for assignment only
    int m_maxshift_high;  //Is used for assignment and assigning in this method
    
    int m_max_ndms;       //This variable is set to 0 in main.cpp and never used until here
    std::vector<std::vector<int>> m_t_processed; //Is allocated in this class, and used elsewhere in the pipeline
    size_t m_t_processed_dim1_size;
    float ***output_buffer; //3D array that contains the output
  };

} //namespace astroaccelerate

#endif /* ASTRO_ACCELERATE_DDTR_STRATEGY */
