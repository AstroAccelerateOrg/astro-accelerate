#ifndef ASTRO_ACCELERATE_AA_ANALYSIS_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_ANALYSIS_STRATEGY_HPP

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <cuda_runtime.h>

#include "aa_strategy.hpp"
#include "aa_analysis_plan.hpp"
#include "aa_params.hpp"
#include "aa_device_MSD_plane_profile.hpp"
#include "aa_device_BC_plan.hpp"
#include "aa_device_info.hpp"

namespace astroaccelerate {

  /** \class aa_analysis_strategy aa_analysis_strategy.hpp "include/aa_analysis_strategy.hpp"
   * \brief Class to configure an analysis strategy.
   * \details An analysis strategy is required for running any pipeline that will run the analysis module aa_compute::module::analysis.
   * \details It is expected behaviour that the configuration values of the strategy may be different than those of the corresponding plan.
   * \author Cees Carels.
   * \date 23 October 2018.
   */  
  class aa_analysis_strategy : public aa_strategy {
  public:
    /** \brief Trivial constructor for aa_analysis_strategy. */
    aa_analysis_strategy() : m_sigma_cutoff(0),
			     m_sigma_constant(0),
			     m_max_boxcar_width_in_sec(0),
			     m_MSD_data_info(0),
			     m_MSD_profile_size_in_bytes(0),
			     m_h_MSD_DIT_width(0),
			     m_candidate_algorithm(aa_analysis_plan::selectable_candidate_algorithm::off),
			     m_enable_sps_baseline_noise(0),
			     m_ready(false) {
      
    }
    
    /**
     * \brief Constructor for aa_analysis_strategy. All parameters must be provided once on construction.
     * \details This constructor is intended to be used when analysis is to be run in isolation, i.e. without running AstroAccelerate's implementation of ddtr.
     * \details In this case, a separate aa_filterbank_metadata must be supplied.
     * \warning At the moment the analysis module does not support this mode, so this constructor still sets m_ready(false) in the initialiser list.
     * \details If the user intends to run the AstroAccelerate implementation of ddtr, then they should use the other non-trivial constructor.
     */
    aa_analysis_strategy(const aa_analysis_plan &analysis_plan, const aa_filterbank_metadata &metadata) : m_metadata(metadata),
													  m_sigma_cutoff(analysis_plan.sigma_cutoff()),
													  m_sigma_constant(analysis_plan.sigma_constant()),
													  m_max_boxcar_width_in_sec(analysis_plan.max_boxcar_width_in_sec()),
													  m_MSD_data_info(0),
													  m_MSD_profile_size_in_bytes(0),
													  m_h_MSD_DIT_width(0),
													  m_candidate_algorithm(analysis_plan.candidate_algorithm()),
													  m_enable_sps_baseline_noise(analysis_plan.enable_sps_baseline_noise()),
													  m_ready(false) {
    }

    /**
     * \brief Constructor for aa_analysis_strategy.
     * \details This constructor is intended to be used when ddtr has also been used.
     * \details Since it uses the aa_filterbank_metadata from ddtr_strategy, the state of aa_analysis_strategy stays consistent with that of aa_ddtr_strategy.
     */
    aa_analysis_strategy(const aa_analysis_plan &analysis_plan) : m_metadata(analysis_plan.ddtr_strategy().metadata()),
								  m_sigma_cutoff(analysis_plan.sigma_cutoff()),
								  m_sigma_constant(analysis_plan.sigma_constant()),
								  m_max_boxcar_width_in_sec(analysis_plan.max_boxcar_width_in_sec()),
								  m_MSD_data_info(0),
								  m_MSD_profile_size_in_bytes(0),
								  m_h_MSD_DIT_width(0),
								  m_candidate_algorithm(analysis_plan.candidate_algorithm()),
								  m_enable_sps_baseline_noise(analysis_plan.enable_sps_baseline_noise()),
								  m_ready(false) {
      if(analysis_plan.ddtr_strategy().configured_for_analysis()) {
	stratagy_MSD(analysis_plan.ddtr_strategy().max_ndms(),
		     analysis_plan.max_boxcar_width_in_sec(),
		     analysis_plan.ddtr_strategy().metadata().tsamp(),
		     analysis_plan.ddtr_strategy().t_processed().at(0).at(0),
		     m_MSD_data_info, m_MSD_profile_size_in_bytes, m_h_MSD_DIT_width);
      }
    }

    /** \returns the name of this strategy type. */
    std::string name() const {
      return "analysis_strategy";
    }
    
    /** \returns the sigma_cutoff determined by the strategy. */
    float sigma_cutoff() const {
      return m_sigma_cutoff;
    }
    
    /** \returns the sigma_constant determined by the strategy. */
    float sigma_constant() const {
      return m_sigma_constant;
    }

    /** \returns the max_boxcar_width_in_sec determined by the strategy. */
    float max_boxcar_width_in_sec() const {
      return m_max_boxcar_width_in_sec;
    }
    
    /** \returns the MSD_data_info determined by the strategy. */ 
    unsigned long int MSD_data_info() const {
      return m_MSD_data_info;
    }
    
    /** \returns the MSD_profile_size_in_bytes determined by the strategy. */
    size_t MSD_profile_size_in_bytes() const {
      return m_MSD_profile_size_in_bytes;
    }
    
    /** \returns the h_MSD_DIT_width determined by the strategy. */
    int h_MSD_DIT_width() const {
      return m_h_MSD_DIT_width;
    }

    /** \returns an integer to indicate whether the candidate algorithm will be enabled or disabled. 0 for off, and 1 for on. */
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

    /** \returns an integer to indicate whether the sps baseline noise reduction algorithm will be enabled or disabled. 0 for off (false), 1 for on (true). */
    int enable_sps_baseline_noise() const {
      return (m_enable_sps_baseline_noise) ? 1 : 0;
    }

    /** \returns The boolean ready state of the strategy, which confirms whether the strategy is valid. */
    bool ready() const {
      return m_ready;
    }
	
    /** \brief Provides a facility to perform any last setup.
     * \details Currently only the ready-state is needed. 
     * \returns A boolean to indicate whether the setup was successful.
     */
    bool setup() {
      return ready();
    }
    
    /** \brief A static member function that prints out the values of the member variables of an instance of aa_analysis_strategy. */
    static bool print_info(const aa_analysis_strategy &strategy) {
      std::cout << "ANALYSIS STRATEGY INFORMATION:" << std::endl;
      std::cout << "analysis sigma_cutoff:\t\t\t" << strategy.sigma_cutoff() <<std::endl;
      std::cout << "analysis sigma_constant:\t\t" << strategy.sigma_constant() << std::endl;
      std::cout << "analysis max_boxcar_width_in_sec:\t" << strategy.max_boxcar_width_in_sec() << std::endl;
      std::cout << "analysis MSD_data_info:\t\t\t" << strategy.MSD_data_info() << std::endl;
      std::cout << "analysis MSD_profile_size_in_bytes:\t" << strategy.MSD_profile_size_in_bytes() << std::endl;
      std::cout << "analysis h_MSD_DIT_width:\t\t" << strategy.h_MSD_DIT_width() << std::endl;
      std::cout << "analysis candidate_algorithm:\t\t" << (strategy.candidate_algorithm() ? "on" : "off") << std::endl;
      std::cout << "analysis sps_baseline_noise:\t\t" << (strategy.enable_sps_baseline_noise() ? "true" : "false") << std::endl;
      return true;
    }
    
  private:
    aa_filterbank_metadata m_metadata; /**< The filterbank metadata for an instance of aa_analysis_strategy. */
    float             m_sigma_cutoff; /**< The sigma_cutoff initially set by aa_analysis_plan, and modified as needed by an instance of aa_analysis_strategy. */
    float             m_sigma_constant; /**< The sigma_constant initially set by the aa_analysis_plan, and modified as needed by an instance of aa_analysis_strategy. */
    float             m_max_boxcar_width_in_sec; /**< Strategy determined m_max_boxcar_width_in_sec. */
    unsigned long int m_MSD_data_info; /**< Strategy determined MSD_data_info. */
    size_t            m_MSD_profile_size_in_bytes; /**< Strategy determined MSD_profile_size_in_bytes. */
    int               m_h_MSD_DIT_width; /**< Strategy determined h_MSD_DIT_width. */
    aa_analysis_plan::selectable_candidate_algorithm m_candidate_algorithm; /**< Flag for selecting candidate algorithm (currently on/off). */
    bool              m_enable_sps_baseline_noise; /**< Flag for enabling/disabling sps_baseline_noise reduction algorithm. */
    bool              m_ready; /**< Ready state for an instance of aa_analysis_strategy. */
    
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

      aa_device_info* mem = aa_device_info::instance();
      size_t free_mem = mem->requested();
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

      // Inform aa_device_memory_manager of the memory that will be required
      // This memory required is given by the allocate_memory_MSD function
      if(!mem->request(maxtimesamples*5.5*sizeof(float))) {
	std::cout << "ERROR:  Could not request memory. " << (unsigned long long)(maxtimesamples*5.5*sizeof(float)) << std::endl;
	return;
      }
      if(!mem->request(sizeof(ushort)*2*maxtimesamples)) {
	std::cout << "ERROR:  Could not request memory. " << (unsigned long long)(sizeof(ushort)*2*maxtimesamples) << std::endl;
	return;
      }
      if(!mem->request(sizeof(float)*MSD_profile_size_in_bytes)) {
	std::cout << "ERROR:  Could not request memory. " << (unsigned long long)(sizeof(float)*MSD_profile_size_in_bytes) << std::endl;
	return;
      }
      
      m_ready = true;
    }
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_ANALYSIS_STRATEGY_HPP
