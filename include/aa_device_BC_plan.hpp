#ifndef ASTRO_ACCELERATE_AA_DEVICE_BC_PLAN_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_BC_PLAN_HPP

namespace astroaccelerate {

  /**
   * \struct PulseDetection_plan aa_device_BC_plan.hpp "include/aa_device_BC_plan.hpp"
   * \brief Data structure for holding pulse detection plan for Box Car (BC) plan.
   * \details The user should not need to interact with this structure directly.
   */

  struct PulseDetection_plan {
    int decimated_timesamples;
    int dtm;
    int iteration;
    int nBoxcars;
    int nBlocks;
    int output_shift;
    int shift;
    int startTaps;
    int unprocessed_samples;
    int total_ut;
  };
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_BC_PLAN_HPP
