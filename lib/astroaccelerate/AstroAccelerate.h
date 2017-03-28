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
 *     The main interface to perform dedispersion and single pulse search
 * @details
 * 	   It fills in a dedispersionPlan (DmTime object) and a candidate list (vector<float>) from the Sps Algorithm
 *
 * @example
 * @code
 *
 * 		// dedispersion strategy
 * 		DedispersionStrategy dedispersion_strategy( list of parameters ...);
 * 		// dedispersed data
 *		DmTime<float> output_buffer(dedispersion_strategy);
 *		// output of sps - assume it's a quarter of the output size
 *		std::vector<float> output_sps;
 *		size_t max_peak_size = (size_t) (dedispersion_strategy.get_nsamp() * dedispersion_strategy.get_nchans()/4);
 *		output_sps.resize(max_peak_size*4);
 *		astroaccelerate::AstroAccelerate<TestParams> astroaccelerate(dedispersion_strategy);
 *		astroaccelerate.run_dedispersion_sps(device_id
 *											,dedispersion_strategy
 *											,input_buffer
 *											,output_buffer
 *											,output_sps
 *											);
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
         * todo: data type for sps and fdas output instead of a vector.
         * sps: dm, timesample, snr, boxcar
         * fdas, f', f, SNR, DM
         * @param device id: the device on which the code should run
         * @param dedispersion_strategy: metadata
         * @param input_buffer: buffer of unsigned short, size number of samples times number of frequency channels
         * @param output_buffer: the dedispersion plan
         * @param output_sps: a vector containing the list of candidates
         * 		  1 candidate is represented by 4 floats: dm, time sample, signal to noise ratio, boxcar filter used
         *
         */
        void run_dedispersion_sps(unsigned device_id
        						  ,DedispersionStrategy &dedispersion_strategy
        						  ,unsigned short *input_buffer
        						  ,DmTime<float> &output_buffer
        						  ,std::vector<float> &output_sps
        						  );
        /*
         * Current state: dd + sps
         * todo: add the following functions once fdas is fully optimized (in progress)
         * void run_dedispersion_fdas();
         * void run_dedispersion_sps_fdas();
         */

        /*
        int get_multi_file() const;
        int get_enable_debug() const;
        int get_enable_analysis() const;
        int get_enable_periodicity() const;
        int get_enable_acceleration() const;
        int get_output_dmt() const;
        int get_enable_zero_dm() const;
        int get_enable_zero_dm_with_outliers() const;
        int get_enable_rfi() const;
        int get_enable_fdas_custom_fft() const;
        int get_enable_fdas_inbin() const;
        int get_enable_fdas_norm() const;
        */
    private:
        /**
         * @brief This function allocates memory for the the gpu arrays based on the dedispersion strategy
         */
        void allocate_memory_gpu(DedispersionStrategy const &);

        /**
         * @brief The number of chunks the data are divided in
         */
        int _num_tchunks;
        /**
         * @brief The number of dm range
         */
        int _range;
        /**
         * @brief The number of frequency channels
         */
        int _nchans;
        /**
         * @brief This value is used to make sure that dms from dm_low to dm_high are used
         */
        int _maxshift;
        /**
         * @brief Time sample value
         */
        float _tsamp;
        /**
         * @brief The maximum number of dm
         */
        int _max_ndms;
        /**
         * @brief The threshold for single pulse detection, multiple of standard deviation
         */
        float _sigma_cutoff;
        /**
         * @brief An array containing the number of dms for each range
         */
        int* _ndms;
        /**
         * @brief An array containing a constant associated with each channel to perform dedispersion algorithm
         */
        float* _dmshifts;
        /**
         * @brief The number of time samples required to search for a dm in each dm range
         */
        int** _t_processed;
        /**
         * @brief An array containing the lowest bound of each dm range
         */
        float* _dm_low;
        /**
         * @brief An array containing the highest bound of each dm range
         */
        float* _dm_high;
        /**
         * @brief An array containing the step size bound of each dm range
         */
        float* _dm_step;
        /**
         * @brief ---
         */
        int* _in_bin;
        /**
         * @brief ---
         */
        int* _out_bin;
        /**
         * @brief Size of the gpu input buffer
         */
        //
        size_t _gpu_input_size;
        /**
         * @brief Gpu input buffer
         */
        unsigned short *_d_input;
        /**
         * @brief Size of the gpu output buffer
         */
        size_t _gpu_output_size;
        /**
         * @brief Gpu output buffer
         */
        float *_d_output;

        // temporary boolean, will disappear eventually
        /**
         * @brief 1 if multiple input files, 0 either
         */
        int _multi_file;
        /**
         * @brief 1 to turn debug on, 0 otherwise
         */
        int _enable_debug;
        /**
         * @brief 1 to turn single pulse search on, 0 otherwise
         */
        int _enable_analysis;
        /**
         * @brief 1 to turn search for periodic objects on, 0 otherwise
         */
        int _enable_periodicity;
        /**
         * @brief 1 to turn fourier domain acceleration search on, 0 otherwise
         */
        int _enable_acceleration;
        /**
         * @brief 1 to turn output dedispspersed data on, 0 otherwise
         */
        int _output_dmt;
        /**
         * @brief 1 to turn analysis with dm=0 on, 0 otherwise
         */
        int _enable_zero_dm;
        /**
         * @brief 1 to turn analysis with dm=0 + outliers on, 0 otherwise
         */
        int _enable_zero_dm_with_outliers;
        /**
         * @brief 1 to turn rfi mitigation on, 0 otherwise
         */
        int _enable_rfi;
        /**
         * @brief 1 to run search with custom FFT, 0 for basic search
         */
        int _enable_fdas_custom_fft;
        /**
         * @brief 1 to perform interbinning on the complex output, 0 otherwise
         */
        int _enable_fdas_inbin;
        /**
         * @brief 1 to perform normalisation, 0 otherwise
         */
        int _enable_fdas_norm;
        /**
         * @brief 1 to use old thresholding code, 0 to use new peak finding code
         */
        int _candidate_algorithm;
};

} // namespace astroaccelerate

#include "detail/AstroAccelerate.cpp"

#endif // ASTROACCELERATE_ASTROACCELERATE_H
