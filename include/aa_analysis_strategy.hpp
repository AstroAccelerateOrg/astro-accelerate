//
//  aa_analysis_strategy.hpp
//  aapipeline
//
//  Created by Cees Carels on Tuesday 23/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef aa_analysis_strategy_hpp
#define aa_analysis_strategy_hpp

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <cuda_runtime.h>

#include "params.hpp"
#include "device_MSD_plane_profile.hpp"
#include "device_BC_plan.hpp"

namespace astroaccelerate {
  
  class aa_analysis_strategy {
  public:
    aa_analysis_strategy() : m_sigma_cutoff(0),
			     m_sigma_constant(0),
			     m_max_boxcar_width_in_sec(0),
			     m_MSD_data_info(0),
			     m_MSD_profile_size_in_bytes(0),
			     m_h_MSD_DIT_width(0),
			     m_candidate_algorithm(aa_analysis_plan::selectable_candidate_algorithm::off),
			     m_enable_sps_baseline_noise(0) {
      
    }
    
    aa_analysis_strategy(const aa_analysis_plan &analysis_plan, const aa_filterbank_metadata &metadata) : m_metadata(metadata),
													  m_sigma_cutoff(analysis_plan.sigma_cutoff()),
													  m_sigma_constant(analysis_plan.sigma_constant()),
													  m_max_boxcar_width_in_sec(analysis_plan.max_boxcar_width_in_sec()),
													  m_MSD_data_info(0),
													  m_MSD_profile_size_in_bytes(0),
													  m_h_MSD_DIT_width(0),
													  m_candidate_algorithm(analysis_plan.candidate_algorithm()),
													  m_enable_sps_baseline_noise(analysis_plan.enable_sps_baseline_noise()) {
      /**
       * Constructor for aa_analysis_strategy.
       * This constructor is intended to be used when analysis is to be run in isolation,
       * that is without having run the AstroAccelerate implementation of ddtr.
       * In this case, a separate aa_filterbank_metadata must be supplied.
       * NOTICE: This functionality is not yet implemented.
       * 
       * If the user intends to run the AstroAccelerate implementation of ddtr, then they should use
       * the constructor that accepts an aa_ddtr_strategy object, as defined below:
       * aa_analysis_strategy::aa_analysis_strategy(const aa_analysis_plan &plan, const aa_ddtr_strategy &ddtr_strategy)
       *
       */
      
    }

    aa_analysis_strategy(const aa_analysis_plan &analysis_plan, const aa_ddtr_strategy &ddtr_strategy) : m_metadata(ddtr_strategy.metadata()),
													 m_sigma_cutoff(analysis_plan.sigma_cutoff()),
													 m_sigma_constant(analysis_plan.sigma_constant()),
													 m_max_boxcar_width_in_sec(analysis_plan.max_boxcar_width_in_sec()),
													 m_MSD_data_info(0),
													 m_MSD_profile_size_in_bytes(0),
													 m_h_MSD_DIT_width(0),
													 m_candidate_algorithm(analysis_plan.candidate_algorithm()),
													 m_enable_sps_baseline_noise(analysis_plan.enable_sps_baseline_noise()) {
      /**
       * Constructor for aa_analysis_strategy.
       * This constructor is intended to be used when ddtr has also been used.
       * Since it uses the aa_filterbank_metadata from ddtr_strategy, the state aa_analysis_strategy
       * stays consistent with that of aa_ddtr_strategy.
       */
      stratagy_MSD(ddtr_strategy.max_ndms(),
		   analysis_plan.max_boxcar_width_in_sec(),
		   ddtr_strategy.metadata().tsamp(),
		   ddtr_strategy.t_processed().at(0).at(0),
		   m_MSD_data_info, m_MSD_profile_size_in_bytes, m_h_MSD_DIT_width);
    }
    
    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }
    
    float sigma_constant() const {
      return m_sigma_constant;
    }

    float max_boxcar_width_in_sec() const {
      return m_max_boxcar_width_in_sec;
    }
    
    unsigned long int MSD_data_info() const {
      return m_MSD_data_info;
    }
    
    size_t MSD_profile_size_in_bytes() const {
      return m_MSD_profile_size_in_bytes;
    }
    
    int h_MSD_DIT_width() const {
      return m_h_MSD_DIT_width;
    }

    int candidate_algorithm() const {
      switch(m_candidate_algorithm) {
      case aa_analysis_plan::selectable_candidate_algorithm::off:
	return 0;
	break;
      case aa_analysis_plan::selectable_candidate_algorithm::on:
	return 1;
	break;
      default:
	return 0;
	break;
      }
      
      return 0;
    }

    int enable_sps_baseline_noise() const {
      return (m_enable_sps_baseline_noise) ? 1 : 0;
    }
    
  private:
    aa_filterbank_metadata m_metadata;
    float             m_sigma_cutoff;
    float             m_sigma_constant;
    float             m_max_boxcar_width_in_sec;
    unsigned long int m_MSD_data_info;
    size_t            m_MSD_profile_size_in_bytes;
    int               m_h_MSD_DIT_width;
    aa_analysis_plan::selectable_candidate_algorithm m_candidate_algorithm;
    bool                m_enable_sps_baseline_noise;
    
    void Create_list_of_boxcar_widths2(std::vector<int> *boxcar_widths, std::vector<int> *BC_widths){
      int DIT_value, DIT_factor, width;
      DIT_value = 1;
      DIT_factor = 2;
      width = 0;
      for(int f=0; f<(int) BC_widths->size(); f++){
	for(int b=0; b<BC_widths->operator[](f); b++){
	  width = width + DIT_value;
	  boxcar_widths->push_back(width);
	}
	DIT_value = DIT_value*DIT_factor;
      }
    }

    void stratagy_MSD(const int &nDMs,
		      const float &max_boxcar_width_in_sec,
		      const float &tsamp,
		      const int &nTimesamples,
		      unsigned long int &maxtimesamples, size_t &MSD_profile_size_in_bytes, int &MSD_nDIT_widths) {

      size_t vals;
      // Calculate the total number of values
      vals = ((size_t) nDMs)*((size_t) nTimesamples);

      size_t free_mem,total_mem;

      cudaMemGetInfo(&free_mem,&total_mem);
      printf("\n----------------------- MSD info ---------------------------\n");
      printf("  Memory required by boxcar filters:%0.3f MB\n",(4.5*vals*sizeof(float) + 2*vals*sizeof(ushort))/(1024.0*1024) );
      printf("  Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
      maxtimesamples = (free_mem*0.95)/((5.5*sizeof(float) + 2*sizeof(ushort)));
      printf("  Max samples: :%lu\n", maxtimesamples);

      int DMs_per_cycle = maxtimesamples/nTimesamples;
      int itemp = (int) (DMs_per_cycle/THR_WARPS_PER_BLOCK);
      DMs_per_cycle = itemp*THR_WARPS_PER_BLOCK;
      printf("\n  DMs_per_cycle: %i",DMs_per_cycle);
	
      int t_BC_widths[10]={PD_MAXTAPS,16,16,16,8,8,8,8,8,8};
      std::vector<int> BC_widths(t_BC_widths,t_BC_widths+sizeof(t_BC_widths)/sizeof(int));
      std::vector<PulseDetection_plan> PD_plan;
      std::vector<int> h_boxcar_widths;
      Create_list_of_boxcar_widths2(&h_boxcar_widths, &BC_widths);

      size_t MSD_DIT_profile_size_in_bytes, workarea_size_in_bytes;
      Get_MSD_plane_profile_memory_requirements(&MSD_profile_size_in_bytes, &MSD_DIT_profile_size_in_bytes, &workarea_size_in_bytes, nTimesamples, nDMs, &h_boxcar_widths);
	
      MSD_nDIT_widths = (int)(max_boxcar_width_in_sec/tsamp);

      printf("\n  Size MSD: %zu \tSize workarea: %f, int: %d\n", MSD_profile_size_in_bytes, workarea_size_in_bytes/1024/1024.0, MSD_nDIT_widths);
      printf("------------------------------------------------------------\n");      
    }
  };

} //namespace astroaccelerate

#endif /* aa_analysis_strategy_hpp */
