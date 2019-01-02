#ifndef ASTRO_ACCELERATE_AA_PIPELINE_RUNNER_HPP
#define ASTRO_ACCELERATE_AA_PIPELINE_RUNNER_HPP

#include <iostream>
#include <vector>

namespace astroaccelerate {

  /**
   * \class aa_pipeline_runner aa_pipeline_runner.hpp "include/aa_pipeline_runner.hpp"
   * \brief Abstract base class for running a pipeline.
   * \details In practice, this class is used for creating base class pointers for pointing to a permitted pipeline instance.
   * \author Cees Carels.
   * \date 30 November 2018.
   */  
  class aa_pipeline_runner {
  public:
    /** \brief Virtual destructor for aa_pipeline_runner. */
    virtual ~aa_pipeline_runner() {
      
    }

    /** \brief Virtual setup method to be implemented by derived classes. */
    virtual bool setup() = 0;
    
    /**
     * \brief Base class virtual methods for running a pipeline.
     * \details In case a derived class does not implement a method, this method will be called.
     */
    virtual bool next() {
      // If a derived class does not implement this method, this method is used.
      std::cout << "ERROR:  The selected operation is not supported on this pipeline." << std::endl;
      return false;
    }

    /** 
     * \brief Base class virtual methods for running a pipeline.
     * \details In case a derived class does not implement a method, this method will be called.
     */
    virtual bool next(std::vector<std::vector<float>> &, int &, std::vector<int> &) {
      // If a derived class does not implement this method, this method is used.
      std::cout << "ERROR:  The selected operation is not supported on this pipeline." << std::endl;
      return false;
    }
  protected:
    
  };
    
} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_PIPELINE_RUNNER_HPP
