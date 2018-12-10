#ifndef ASTRO_ACCELERATE_AA_PIPELINE_RUNNER_HPP
#define ASTRO_ACCELERATE_AA_PIPELINE_RUNNER_HPP

#include <iostream>
#include <vector>

namespace astroaccelerate {
  
  class aa_pipeline_runner {
  public:
    virtual ~aa_pipeline_runner() {
      
    }

    virtual bool setup() = 0;
    
    virtual bool next() {
      // If a derived class does not implement this method, this method is used.
      std::cout << "ERROR:  The selected operation is not supported on this pipeline." << std::endl;
      return false;
    }
    virtual bool next(std::vector<float> &, int &, std::vector<int> &) {
      // If a derived class does not implement this method, this method is used.
      std::cout << "ERROR:  The selected operation is not supported on this pipeline." << std::endl;
      return false;
    }
  protected:
    
  };
    
} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_PIPELINE_RUNNER_HPP
