#ifndef SKA_ASTROACCELERATE_SPS_SPS_H
#define SKA_ASTROACCELERATE_SPS_SPS_H

#include "DedispersionPlan.h"
#include "IOData.h"
#include "UserInput.h"

#include "../AstroAccelerate/device_init.h"
#include "../AstroAccelerate/device_load_data.h"
#include "../AstroAccelerate/device_zero_dm.h"
#include "../AstroAccelerate/device_corner_turn.h"
#include "../AstroAccelerate/device_dedisperse.h"
#include "../AstroAccelerate/device_bin.h"
#include "../AstroAccelerate/device_save_data.h"
#include "../AstroAccelerate/host_write_file.h"
#include "../AstroAccelerate/host_analysis.h"

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
        void operator()(unsigned device_id, IOData &, DedispersionPlan &,
                        UserInput const &, SpsHandler, DmHandler);

    private:
};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
#include "detail/Sps.cpp"

#endif // SKA_ASTROACCELERATE_SPS_SPS_H
