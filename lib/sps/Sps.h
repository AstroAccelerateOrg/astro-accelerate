#ifndef SKA_ASTROACCELERATE_SPS_SPS_H
#define SKA_ASTROACCELERATE_SPS_SPS_H

#include "DedispersionPlan.h"
#include "IOData.h"
#include "UserInput.h"


namespace ska {
namespace astroaccelerate {
namespace sps {

/**
 * @brief
 *     The main interface to perform single pulse search and dedispersion
 * @details
 * @tparam SpsParameterType the optimisation compile time parameters
 *         for the Sps algorithm.....
 *         nthreads...
 *
 * @example
 * @code
 *
 *
 * @endcode
 * @endexample
 */

template<typename SpsParameterType>
class Sps
{
    public:
        Sps();
        ~Sps();

        /**
         * @brief perform dedispersion and an sps search
         */
        void operator()(unsigned device_id, IOData &, DedispersionPlan &, UserInput const &);

    private:

};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
#include "detail/Sps.cpp"

#endif // SKA_ASTROACCELERATE_SPS_SPS_H
