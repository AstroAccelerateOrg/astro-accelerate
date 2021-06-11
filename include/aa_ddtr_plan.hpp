#ifndef ASTRO_ACCELERATE_AA_DDTR_PLAN_HPP
#define ASTRO_ACCELERATE_AA_DDTR_PLAN_HPP

#include <stdio.h>
#include <vector>

namespace astroaccelerate {

  /**
   * \class aa_ddtr_plan aa_ddtr_plan.hpp "include/aa_ddtr_plan.hpp"
   * \brief Class to set a ddtr plan.
   * \details A ddtr plan is required in order to create a ddtr strategy.
   * \author Cees Carels.
   * \date 23 October 2018.
   */

  class aa_ddtr_plan {
  public:
    /** \brief Constructor for aa_ddtr_plan. m_power must be default initialised to 2.0. */
    aa_ddtr_plan() : m_power(2.0), m_enable_msd_baseline_noise(false) {
      
    }
  
    /** \struct dm
     * \brief A struct for specifying dm ranges and binning.
     * \deprecated dm.outBin is deprecated.
     */
    struct dm {
      float low;   /**< The low value of a dm range. */
      float high;  /**< The high value of a dm range. */
      float step;  /**< The step size of a dm range. */
      int   inBin; /**< The binning factor of a dm range. */
      int   outBin;/**< \deprecated Variable is deprecated. */
    };

    /** \brief Add a single dm range by setting all required values from struct dm individually.
     * \return True if the operation is successful, false otherwise. 
     */
    bool add_dm(const float &low, const float &high, const float &step, const int &inBin, const int &outBin) {
      const dm tmp = {low, high, step, inBin, outBin};
      m_user_dm.push_back(std::move(tmp));
      return true;
    }

    /** \brief Add a single dm range by providing a struct dm object directly.
     * \return True if the operation is successful, false otherwise.
     */
    bool add_dm(aa_ddtr_plan::dm &DM) {
      m_user_dm.push_back(std::move(DM));
      return true;
    }

    /** \brief Add a std::vector of struct dm objects directly.
     * \return True if the operation is successful, false otherwise. 
     */
    bool add_dm(std::vector<aa_ddtr_plan::dm> &DM) {
      m_user_dm.insert(end(m_user_dm), begin(DM), end(DM));
      return true;
    }

    /** \brief Get the total number of dm ranges that were supplied.
     * \return The number of supplied dm ranges.
     */
    size_t range() const {
      return m_user_dm.size();
    }

    /**
     * \brief Get the supplied dm range at index i.
     * \details Using add_dm, the supplied dm ranges can be accessed.
     * \return The dm object at index i.
     * \warning Attemping to access out of bounds will return a dm object with all values set to 0.
     */
    const dm user_dm(const size_t &i) const {
      if(i < m_user_dm.size()) {
	return m_user_dm.at(i);
      }
      else {
	aa_ddtr_plan::dm empty_dm = {0.0, 0.0, 0.0, 0, 0};
	return empty_dm;
      }
    }

    /**
     * \brief Set the ddtr power, which is needed to calculate the aa_ddtr_strategy instance needed for dedispersion.
     */
    bool set_power(const float &power) {
      m_power = power;
      return true;
    }

    /** \returns The user set power. */
    float power() const {
      return m_power;
    }

    /**
     * \brief Set flag to enable or disable msd_baseline_noise reduction algorithm.
     * \details At present, this flag has no effect.
     */
    bool set_enable_msd_baseline_noise(const bool flag) {
      m_enable_msd_baseline_noise = flag;
      return true;
    }

    /**
     * \returns The flag that enables or disables msd_baseline_noise reduction.
     * \details At the moment, this setting has no effect.
     */
    bool enable_msd_baseline_noise() const {
      return m_enable_msd_baseline_noise;
    }
	
    /**
     * \brief Set the custom bandpass normalization values for zerodm filtering.
     */
    bool bind_bandpass_normalization(const float *custom_bandpass_normalization, int bandpass_size) {
      if(bandpass_size>0){
        bandpass_normalization.resize(bandpass_size);
        std::copy( custom_bandpass_normalization, custom_bandpass_normalization + bandpass_size, bandpass_normalization.begin() );
		return true;
	  }
	  else {
        return false;
	  }
    }
	
    /**
     * \brief Set the custom bandpass normalization values for zerodm filtering.
     */
	bool bind_bandpass_normalization_vector(std::vector<float> &custom_bandpass_normalization) {
      bandpass_normalization = custom_bandpass_normalization;
      return true;
	}
	
    /**
     * \returns the size of the custom bandpass normalization array.
	 */
	size_t bandpass_normalization_size() const {
		return(bandpass_normalization.size());
	}
	
    /**
     * \returns the pointer to the custom bandpass normalization array.
	 */
	const float* bandpass_data_pointer() const {
		if(bandpass_normalization.size()>0){
			return(bandpass_normalization.data());
		}
		else return NULL;
	}
    
  private:
    std::vector<dm> m_user_dm; /**< Storage for all supplied dm properties. */
    float m_power;
    bool m_enable_msd_baseline_noise; /** Flag to enable or disable msd_baseline_noise reduction algorithm. */
	std::vector<float> bandpass_normalization;
  
  };
} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_DDTR_PLAN_HPP
