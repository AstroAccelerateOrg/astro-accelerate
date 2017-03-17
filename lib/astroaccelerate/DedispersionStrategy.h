#ifndef ASTROACCELERATE_SPS_DEDISPERSIONSTRATEGY_H
#define ASTROACCELERATE_SPS_DEDISPERSIONSTRATEGY_H

#include "../headers/params.h"
#include "../headers/headers_mains.h"
#include "../headers/host_help.h"

#include <stdio.h>

namespace astroaccelerate {

    /**
     * @brief  Dedispersion Strategy
     *
     * @details This object carries the dedispersion strategy
     *
     */
class DedispersionStrategy
{
    public:
        /**
        *  @brief Default constructor
        */
        DedispersionStrategy();

        /**
        *  @brief Destructor
        */
        ~DedispersionStrategy();

        /**
        *  @brief Getters
        */
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
        int get_nboots() const ;
        int get_ntrial_bins() const;
        int get_navdms() const;
        float get_narrow() const;
        float get_aggression() const;
        int get_nsearch() const;
        float get_power() const;
        float get_sigma_cutoff() const;
        float get_sigma_constant() const;
        float get_max_boxcar_width_in_sec() const;
        float get_wide() const;
        int get_range() const;
        float* get_user_dm_low() const;
        float* get_user_dm_high() const;
        float* get_user_dm_step() const;
        int* get_in_bin() const;
        int* get_out_bin() const;
        //
        int get_maxshift() const;
        float* get_dm_low() const;
        float* get_dm_high() const;
        float* get_dm_step() const;
        float* get_dmshifts() const;
        int* get_ndms() const ;
        int get_max_ndms() const;
        int get_total_ndms() const;
        float get_max_dm() const;
        int** get_t_processed() const;
        int get_nbits() const;
        int get_nifs() const;
        float get_tstart() const;
        float get_tsamp() const;
        int get_nsamp() const;
        int get_nsamples() const;
        int get_max_samps() const;
        int get_nchans() const;
        float get_fch1() const;
        float get_foff() const;
        unsigned int get_num_tchunks() const;

        /**
        *  @brief   Get the user input
        *  @details This function read the user input file and stores it
        *
        */
        void get_user_input(FILE** fp
        					,int argc
        					,char *argv[]);
        /**
         * @brief Read telescope parameters from the header of the input file
         */
        void get_file_data(FILE **);

        /**
         * @brief Computes the dedispersion strategy
         */
        void make_strategy(size_t);

    private:
        // user input
        int _multi_file;
        int _enable_debug;
        int _enable_analysis;
        int _enable_periodicity;
        int _enable_acceleration;
        int _output_dmt;
        int _enable_zero_dm;
        int _enable_zero_dm_with_outliers;
        int _enable_rfi;
        int _enable_fdas_custom_fft;
        int _enable_fdas_inbin;
        int _enable_fdas_norm;
        int _nboots;
        int _ntrial_bins;
        int _navdms;
        float _narrow;
        float _aggression;
        int _nsearch;
        float _power;
        float _sigma_cutoff;
        float _sigma_constant;
        float _max_boxcar_width_in_sec;
        float _wide;
        int _range;
        float* _user_dm_low;
        float* _user_dm_high;
        float* _user_dm_step;
        int* _in_bin;
        int* _out_bin;
        // dedispersion strategy
        int _maxshift;
        float* _dm_low;
        float* _dm_high;
        float* _dm_step;
        float* _dmshifts;
        int* _ndms;
        int _max_ndms;
        int _total_ndms;
        float _max_dm;
        int** _t_processed;
        int _nbits;
        int _nifs;
        float _tstart;
        float _tsamp;
        int _nsamp;
        int _nsamples;
        int _max_samps;
        int _nchans;
        float _fch1;
        float  _foff;
        unsigned int _num_tchunks;
};

} // namespace astroaccelerate


#endif // ASTROACCELERATE_DEDISPERSIONSTRATEGY_H
