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
        Sps(IOData &,
        	DedispersionPlan &,
        	UserInput &);
        ~Sps();

        /**
         * @brief perform dedispersion and an sps search
         */
        void operator()(unsigned device_id, IOData &, DedispersionPlan &, UserInput &);

    private:
        int _num_tchunks;
        int _range;
        int _nchans;
        int _maxshift;
        float _tsamp;
        int _max_ndms;
        float _sigma_cutoff;
        size_t _gpu_outputsize;
        int* _ndms;
        float* _dmshifts;
        int** _t_processed;
        float* _dm_low;
        float* _dm_high;
        float* _dm_step;
        int* _in_bin;
        int* _out_bin;
        unsigned short* _input_buffer;
        float *** _output_buffer;
		unsigned short* _d_input;
		float* _d_output;

};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
#include "detail/Sps.cpp"

#endif // SKA_ASTROACCELERATE_SPS_SPS_H
