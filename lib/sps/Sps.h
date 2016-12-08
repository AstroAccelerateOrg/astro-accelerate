#ifndef SKA_ASTROACCELERATE_SPS_SPS_H
#define SKA_ASTROACCELERATE_SPS_SPS_H

#include "DedispersionPlan.h"
#include "IOData.h"
#include "../AstroAccelerate/device_init.h"

namespace ska {
namespace astroaccelerate {
namespace sps {

/**
 * @brief
 *     The main interface to perform single pulse search and dedispersion
 * @details
 * @tparam SpsParameterType the optimisiation compile time parameters
 *         for the Sps algortihm.....
 *         nthreads...
 *
 * @example
 * @code
 *
 * Sps<MyParams> sps;
 * sps(in, out, );
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
        template<typename SpsHandler, typename DmHandler>
        void operator()(unsigned device_id, IOData &,
                        DedispersionPlan &, SpsHandler, DmHandler);

    private:
};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
#include "detail/Sps.cpp"

#endif // SKA_ASTROACCELERATE_SPS_SPS_H
