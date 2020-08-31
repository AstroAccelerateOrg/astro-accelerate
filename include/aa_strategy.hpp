#ifndef ASTRO_ACCELERATE_AA_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_STRATEGY_HPP

#include <stdio.h>
#include <string>

namespace astroaccelerate {

  /** 
   * \class aa_strategy aa_strategy.hpp "include/aa_strategy.hpp"
   * \brief Abstract base class for all strategy classes.
   * \details All strategy classes must adhere to this interface so that the API can setup and ready all component.
   * \author Cees Carels.
   * \date 1 November 2018.
   */
  class aa_strategy {
  public:
    virtual bool        setup()       = 0; /** \brief Setup the strategy. \returns A boolean flag to indicate whether the operation was successful. */
    virtual std::string name()  const = 0; /** \returns The name of the component. */
  private:
  
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_STRATEGY_HPP
