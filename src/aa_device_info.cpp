#include "aa_device_info.hpp"

namespace astroaccelerate {

  aa_device_info *aa_device_info::m_instance = 0;
  bool aa_device_info::is_init = false;
  std::vector<aa_device_info::aa_card_info> aa_device_info::m_card_info;
  
} // namespace astroaccelerate
