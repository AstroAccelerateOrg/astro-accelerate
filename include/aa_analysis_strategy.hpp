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
   * \details An analysis strategy is required for running any pipeline that will run the analysis component aa_pipeline::component::analysis.
   * \details It is expected behaviour that the configuration values of the strategy may be different than those of the corresponding plan.
   * \author AstroAccelerate
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
				 m_selected_device(NULL),
			     m_candidate_algorithm(aa_analysis_plan::selectable_candidate_algorithm::peak_find),
			     m_enable_msd_baseline_noise(0),
			     m_ready(false) {
      
    }
    
    /**
     * \brief Constructor for aa_analysis_strategy. All parameters must be provided once on construction.
     * \details This constructor is intended to be used when analysis is to be run in isolation, i.e. without running AstroAccelerate's implementation of ddtr.
     * \details In this case, a separate aa_filterbank_metadata must be supplied.
     * \warning At the moment the analysis component does not support this mode, so this constructor still sets m_ready(false) in the initialiser list.
     * \details If the user intends to run the AstroAccelerate implementation of ddtr, then they should use the other non-trivial constructor.
     */
    aa_analysis_strategy(const aa_analysis_plan &analysis_plan, const aa_filterbank_metadata &metadata, aa_device_info *selected_device) : m_metadata(metadata),
													  m_sigma_cutoff(analysis_plan.sigma_cutoff()),
													  m_sigma_constant(analysis_plan.sigma_constant()),
													  m_max_boxcar_width_in_sec(analysis_plan.max_boxcar_width_in_sec()),
													  m_MSD_data_info(0),
													  m_MSD_profile_size_in_bytes(0),
													  m_h_MSD_DIT_width(0),
													  m_selected_device(selected_device),
													  m_candidate_algorithm(analysis_plan.candidate_algorithm()),
													  m_enable_msd_baseline_noise(analysis_plan.enable_msd_baseline_noise()),
													  m_ready(false) {
    }

    /**
     * \brief Constructor for aa_analysis_strategy.
     * \details This constructor is intended to be used when ddtr has also been used.
     * \details Since it uses the aa_filterbank_metadata from ddtr_strategy, the state of aa_analysis_strategy stays consistent with that of aa_ddtr_strategy.
     */
    aa_analysis_strategy(const aa_analysis_plan &analysis_plan, aa_device_info *selected_device) : m_metadata(analysis_plan.ddtr_strategy().metadata()),
								  m_sigma_cutoff(analysis_plan.sigma_cutoff()),
								  m_sigma_constant(analysis_plan.sigma_constant()),
								  m_max_boxcar_width_in_sec(analysis_plan.max_boxcar_width_in_sec()),
								  m_MSD_data_info(0),
								  m_MSD_profile_size_in_bytes(0),
								  m_h_MSD_DIT_width(0),
								  m_selected_device(selected_device),
								  m_candidate_algorithm(analysis_plan.candidate_algorithm()),
								  m_enable_msd_baseline_noise(analysis_plan.enable_msd_baseline_noise()),
								  m_ready(false) {
      
      bool ready_to_configure = false;
      if((m_sigma_cutoff > 0)
	 &&
	 (m_sigma_constant > 0)
	 &&
	 (m_max_boxcar_width_in_sec > 0)) {
	ready_to_configure = true;
      }
      
      if((ready_to_configure) && (analysis_plan.ddtr_strategy().configured_for_analysis())) {
	stratagy_MSD(
		analysis_plan.ddtr_strategy().max_ndms(),
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
			case aa_analysis_plan::selectable_candidate_algorithm::peak_find:
				return 0;
				break;
			case aa_analysis_plan::selectable_candidate_algorithm::threshold:
				return 1;
				break;
			case aa_analysis_plan::selectable_candidate_algorithm::peak_filtering:
				return 2;
				break;
			default:
				return 0;
				break;
		}
		return 0;
    }

    /** \returns an integer to indicate whether the msd baseline noise reduction algorithm will be enabled or disabled. 0 for off (false), 1 for on (true). */
    int enable_msd_baseline_noise() const {
      return (m_enable_msd_baseline_noise) ? 1 : 0;
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
      LOG(log_level::dev_debug, "ANALYSIS STRATEGY INFORMATION:");
      LOG(log_level::dev_debug, "analysis sigma_cutoff:\t\t\t" + std::to_string(strategy.sigma_cutoff()));
      LOG(log_level::dev_debug, "analysis sigma_constant:\t\t\t" + std::to_string(strategy.sigma_constant()));
      LOG(log_level::dev_debug, "analysis max_boxcar_width_in_sec:\t" + std::to_string(strategy.max_boxcar_width_in_sec()));
      LOG(log_level::dev_debug, "analysis MSD_data_info:\t\t\t" + std::to_string(strategy.MSD_data_info()));
      LOG(log_level::dev_debug, "analysis MSD_profile_size_in_bytes:\t" + std::to_string(strategy.MSD_profile_size_in_bytes()));
      LOG(log_level::dev_debug, "analysis h_MSD_DIT_width:\t\t" + std::to_string(strategy.h_MSD_DIT_width()));
      LOG(log_level::dev_debug, "analysis candidate_algorithm:\t\t" + (strategy.candidate_algorithm() ? std::string("on") : std::string("off")));
      LOG(log_level::dev_debug, "analysis msd_baseline_noise:\t\t" + (strategy.enable_msd_baseline_noise() ? std::string("true") : std::string("false")));
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
	aa_device_info    *m_selected_device;
    aa_analysis_plan::selectable_candidate_algorithm m_candidate_algorithm; /**< Flag for selecting candidate algorithm (currently on/off). */
    bool              m_enable_msd_baseline_noise; /**< Flag for enabling/disabling msd_baseline_noise reduction algorithm. */
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

    void stratagy_MSD(
			const int &nDMs,
			const int &nTimesamples,
			unsigned long int &maxtimesamples, size_t &MSD_profile_size_in_bytes, int &MSD_nDIT_widths) {

      size_t vals;
      // Calculate the total number of values
      vals = ((size_t) nDMs)*((size_t) nTimesamples);

      //aa_device_info& mem = aa_device_info::instance();
      size_t free_mem = m_selected_device->free_memory();
      printf("\n----------------------- MSD info ---------------------------\n");
      printf("  Memory required by single pulse detection:%0.3f MB\n",(5.5*vals*sizeof(float) + 2*vals*sizeof(ushort))/(1024.0*1024) );
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
	
      MSD_nDIT_widths = h_boxcar_widths.size();

      printf("  Number of boxcar filter steps: %zu;\n", h_boxcar_widths.size());
      printf("  Memory required for MSD results for all boxcar filters: %zu bytes;\n", MSD_profile_size_in_bytes);
      printf("  Size of MSD temporary workarea: %0.3f MB;\n", workarea_size_in_bytes/1024/1024.0);
      printf("------------------------------------------------------------\n");

      // Inform aa_device_memory_manager of the memory that will be required
      // This memory required is given by the allocate_memory_MSD function
      if(!m_selected_device->request_memory(maxtimesamples*5.5*sizeof(float))) {
        std::cout << "ERROR:  Could not request memory. " << (unsigned long long)(maxtimesamples*5.5*sizeof(float)) << std::endl;
        return;
      }
      if(!m_selected_device->request_memory(sizeof(ushort)*2*maxtimesamples)) {
        std::cout << "ERROR:  Could not request memory. " << (unsigned long long)(sizeof(ushort)*2*maxtimesamples) << std::endl;
        return;
      }
      if(!m_selected_device->request_memory(sizeof(float)*MSD_profile_size_in_bytes)) {
        std::cout << "ERROR:  Could not request memory. " << (unsigned long long)(sizeof(float)*MSD_profile_size_in_bytes) << std::endl;
        return;
      }
	  
      
      m_ready = true;
    }
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_ANALYSIS_STRATEGY_HPP
