#ifndef ASTRO_ACCELERATE_AA_DEVICE_PERIODS_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_PERIODS_HPP

#include "aa_periodicity_strategy.hpp"
#include "aa_periodicity_candidates.hpp"

namespace astroaccelerate {

extern void GPU_periodicity(aa_periodicity_strategy &PSR_strategy, float ***output_buffer, aa_periodicity_candidates &Power_Candidates, aa_periodicity_candidates &Interbin_Candidates);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_PERIODS_HPP








