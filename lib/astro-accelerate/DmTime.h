#ifndef ASTROACCELERATE_CUDA_DMTIME_H
#define ASTROACCELERATE_CUDA_DMTIME_H

#include "DedispersionStrategy.h"
#include <vector>
#include <assert.h>

namespace astroaccelerate {

/**
 * @brief
 *    Class to encapsulate DmTime
 *    Credit goes to Chris Williams
 * @details
 *
 */

template<typename ValueType>
class DmTime
{
    public:
        DmTime(DedispersionStrategy const&);
        ~DmTime();
        std::size_t number_of_dm_ranges() const;
        std::size_t output_size() const;
        std::vector<std::size_t> const& nsamples() const;
        float** operator[](std::size_t dm_range) { return _data[dm_range]; };

    private:
        /**
         * @brief Size of the dedispersion plan
         */
        std::size_t _ouput_size;

        /**
         * @brief The number of time samples for the set of DM trials in each range
         */
        std::vector<std::size_t> _nsamples;
        /**
         * @brief Number of dms for each dm range
         */
        std::vector<std::size_t> _ndms;
        /**
         * @brief Number of dm range
         */
        std::size_t _range;
        /**
         * @brief Dedispersion plan
         * @details time * dms * power of signal
         */
        ValueType*** _data; // contiguous memory block
};

}

#include "detail/DmTime.cpp"
#endif // ASTROACCELERATE_CUDA_DMTIME_H
