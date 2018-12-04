//
//  aa_strategy.hpp
//  aapipeline
//
//  Created by Cees Carels on Thursday 01/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_STRATEGY_HPP
#define ASTRO_ACCELERATE_STRATEGY_HPP

#include <stdio.h>

class aa_strategy {
public:
  virtual bool        setup()       = 0;
  virtual std::string name()  const = 0;
private:
  
};

#endif /* ASTRO_ACCELERATE_STRATEGY_HPP */
