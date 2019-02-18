#include "aa_welcome_notice.hpp"

#include "aa_log.hpp"
namespace astroaccelerate {
  void welcome_notice() {
    LOG(log_level::notice, "=== STARTING ASTRO-ACCELERATE ===");
    LOG(log_level::notice, "Please cite the references provided in the repository.");
    LOG(log_level::notice, "https://github.com/AstroAccelerateOrg/astro-accelerate#publications");
  }
} // namespace astroaccelerate
