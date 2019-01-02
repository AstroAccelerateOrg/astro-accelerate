#ifndef ASTRO_ACCELERATE_AA_DDTR_PLAN_HPP
#define ASTRO_ACCELERATE_AA_DDTR_PLAN_HPP

#include <stdio.h>
#include <vector>

/**
 * \class aa_ddtr_plan aa_ddtr_plan.hpp "include/aa_ddtr_plan.hpp"
 * \brief Class to set a ddtr plan.
 * \details A ddtr plan is required in order to create a ddtr strategy.
 * \author Cees Carels.
 * \date 23 October 2018.
 */

class aa_ddtr_plan {
public:
  /** \brief Constructor for aa_ddtr_plan. */
  aa_ddtr_plan() {
      
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
    
private:
  std::vector<dm> m_user_dm; /**< Storage for all supplied dm properties. */
  
};

#endif /* ASTRO_ACCELERATE_AA_DDTR_PLAN_HPP */
