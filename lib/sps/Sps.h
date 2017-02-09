#ifndef SKA_ASTROACCELERATE_SPS_SPS_H
#define SKA_ASTROACCELERATE_SPS_SPS_H

#include "DedispersionPlan.h"
#include "InputData.h"
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
        Sps(InputData &,
            DedispersionPlan &,
            UserInput &);
        ~Sps();

        /**
         * @brief allocate memory for the cpu output and the gpu
         */
        void allocate_memory_cpu_output(DedispersionPlan const &);
        void allocate_memory_gpu(DedispersionPlan const &);

        /**
         * @brief perform dedispersion and an sps search
         */
        void operator()(unsigned device_id, InputData &, DedispersionPlan &, UserInput &, size_t gpu_memory);

    private:

        int _num_tchunks;
        int _range;
        int _nchans;
        int _maxshift;
        float _tsamp;
        int _max_ndms;
        float _sigma_cutoff;
        int* _ndms;
        float* _dmshifts;
        int** _t_processed;
        float* _dm_low;
        float* _dm_high;
        float* _dm_step;
        int* _in_bin;
        int* _out_bin;
        unsigned short* _input_buffer;
        //
        size_t  _output_size;
        float *** _output_buffer;
        size_t          _gpu_input_size;
        unsigned short  *_d_input;
        size_t  _gpu_output_size;
        float   *_d_output;

};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
#include "detail/Sps.cpp"

#endif // SKA_ASTROACCELERATE_SPS_SPS_H
