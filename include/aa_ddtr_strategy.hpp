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

class aa_ddtr_strategy : public aa_strategy {
public:
    aa_ddtr_strategy();
    aa_ddtr_strategy(const aa_ddtr_plan &plan, const aa_filterbank_metadata &metadata, const size_t &free_memory, const bool &enable_analysis);
    bool setup();
    
    int maxshift() const {
        return m_maxshift;
    }
    
    int** t_processed() const {
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
    
private:
    bool strategy(const aa_ddtr_plan &plan, const aa_filterbank_metadata &metadata, const size_t &free_memory, const bool &enable_analysis);
    bool m_ready;
    void allocate_memory_cpu_output();
    
    bool is_setup;  //Has setup been called already?
    
    std::vector<int> m_ndms;
    std::vector<float> m_dmshifts;
    std::vector<aa_ddtr_plan::dm> str_dm;
    int m_maxshift;       //Is used for assignment and assigning
    int m_num_tchunks;    //Is used for assignment only
    int m_total_ndms;     //Is used for assignment only
    float m_max_dm;       //Is used for assignment only
    int m_maxshift_high;  //Is used for assignment and assigning in this method
    
    int m_max_ndms;       //This variable is set to 0 in main.cpp and never used until here
    int **m_t_processed; //Is allocated in this class, and used elsewhere in the pipeline
    float ***output_buffer; //3D array that contains the output
};

#endif /* ASTRO_ACCELERATE_DDTR_STRATEGY */
