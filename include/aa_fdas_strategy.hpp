//
//  aa_fdas_strategy.hpp
//  aapipeline
//
//  Created by Cees Carels on Monday 03/12/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP

#include <stdio.h>

#include "aa_strategy.hpp"
#include "aa_fdas_plan.hpp"

namespace astroaccelerate {

  /**
   * Class for aa_fdas_strategy, used to configure the fourier domain accelerated search (fdas).
   */

  class aa_fdas_strategy : public aa_strategy {
  public:
    aa_fdas_strategy() : m_ready(false) {
      
    }
    
    aa_fdas_strategy(const aa_fdas_plan &fdas_plan) : m_ready(false) {
      
    }

    bool setup() {
      return false;
    }

    bool ready() const {
      return m_ready;
    }
  private:
    bool m_ready;

  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_FDAS_STRATEGY_HPP
