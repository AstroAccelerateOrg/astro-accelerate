#ifndef ASTROACCELERATE_ASTROACCELERATE_H
#define ASTROACCELERATE_ASTROACCELERATE_H

#include "DedispersionStrategy.h"
#include "DmTime.h"
//
#include "../headers/headers_mains.h"
#include "../headers/device_bin.h"
#include "../headers/device_init.h"
#include "../headers/device_dedisperse.h"
#include "../headers/device_dedispersion_kernel.h"
#include "../headers/device_zero_dm.h"
#include "../headers/device_zero_dm_outliers.h"
#include "../headers/device_rfi.h"
// sps
#include "../headers/device_BLN.h" //Added by KA
#include "../headers/device_SPS_inplace_kernel.h" //Added by KA
#include "../headers/device_SPS_inplace.h" //Added by KA
#include "../headers/device_MSD_grid.h" //Added by KA
#include "../headers/device_MSD_plane.h" //Added by KA
#include "../headers/device_MSD_limited.h" //Added by KA
#include "../headers/device_SNR_limited.h" //Added by KA
#include "../headers/device_threshold.h" //Added by KA
#include "../headers/device_single_FIR.h" //Added by KA
#include "../headers/device_analysis.h" //Added by KA
#include "../headers/device_peak_find.h" //Added by KA
//
#include "../headers/device_load_data.h"
#include "../headers/device_corner_turn.h"
#include "../headers/device_save_data.h"
#include "../headers/host_acceleration.h"
#include "../headers/host_allocate_memory.h"
#include "../headers/host_analysis.h"
#include "../headers/host_periods.h"
#include "../headers/host_debug.h"
#include "../headers/host_get_file_data.h"
#include "../headers/host_get_recorded_data.h"
#include "../headers/host_get_user_input.h"
#include "../headers/host_help.h"
#include "../headers/host_rfi.h"
#include "../headers/host_stratagy.h"
#include "../headers/host_write_file.h"
// fdas
#include "../headers/device_acceleration_fdas.h"
#include "../headers/host_main_function.h"
//
#include "../headers/params.h"
#include "timer.h"

#include <vector>

namespace astroaccelerate {

/**
 * @brief
 *     The main interface to perform dedispersion, single pulse search and fdas
 * @details
 *
 * @example
 * @code
 *
 *
 * @endcode
 * @endexample
 */

template<typename AstroAccelerateParameterType>
class AstroAccelerate
{
	public:
        AstroAccelerate(DedispersionStrategy &);
        ~AstroAccelerate();

        /**
         * @brief perform dedispersion , sps search, fdas
         * todo: object for sps and fdas output instead of a vector.
         * sps: dm, timesample, snr, boxcar
         * fdas, f', f, SNR, DM
         */
        void run_dedispersion_sps_fdas(unsigned device_id
        						  ,DedispersionStrategy &
        						  //,size_t gpu_memory
        						  ,unsigned short *input_buffer
        						  ,DmTime<float> &output_buffer
        						  ,std::vector<float> &output_sps
        						  );
        /*
         * todo: define once got rid of fdas linking problem ...
         */
/*        void run_dedispersion_sps();
          void run_dedispersion_fdas();
          void run_dedispersion_sps_fdas();
*/
    private:
        /**
         * @brief allocate memory for the the gpu arrays
         */
        void allocate_memory_gpu(DedispersionStrategy const &);

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
        //
        size_t _gpu_input_size;
        unsigned short *_d_input;
        size_t _gpu_output_size;
        float *_d_output;

};

} // namespace astroaccelerate

#include "detail/AstroAccelerate.cpp"

#endif // ASTROACCELERATE_ASTROACCELERATE_H
