#ifndef ASTROACCELERATE_SPS_DEDISPERSIONSTRATEGY_H
#define ASTROACCELERATE_SPS_DEDISPERSIONSTRATEGY_H

#include "../headers/params.h"
#include "../headers/headers_mains.h"
#include "../headers/host_help.h"

#include <stdio.h>
#include <vector>

namespace astroaccelerate {

    /**
     * @brief  Dedispersion Strategy
     *
     * @details This object carries the dedispersion strategy
     *
     */
class DedispersionStrategy
{
	friend class DedispersionStrategyFile;
    public:
        /**
        *  @brief Default constructor
        */
        DedispersionStrategy();
        /**
         *  @brief Parameterized constructor
         */
        DedispersionStrategy(float* const user_dm_low
        			,float* const user_dm_high
        			,float* const user_dm_step
        			,int* const in_bin
        			,int* const out_bin
        			,size_t gpu_memory
        			,int power
        			,int range
        			,int nchans
        			,int nsamples
        			,int nsamp
        			,int nifs
        			,int nbits
        			,float tsamp
        			,float tstart
        			,float sigma_cutoff
				,float sigma_constant
				,float max_boxcar_width_in_sec
				,float narrow
				,float wide
				,int nboots
				,int navdms
				,int ntrial_bins
				,int nsearch
				,float aggression
                                ,std::vector<float> const & bin_frequencies
				);
        /**
        *  @brief Destructor
        */
        ~DedispersionStrategy();

        /**
        *  @brief Getters
        */
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
        unsigned int get_num_tchunks() const;

        void resize(size_t number_of_samples, size_t gpu_memory);

    private:
        /**
         * @brief Computes the dedispersion strategy
         *
         */
        void make_strategy(size_t gpu_memory);
        
 
	// user input
        /**
         * @brief ---
         */
        int _nboots;
        /**
         * @brief ---
         */
        int _ntrial_bins;
        /**
         * @brief ---
         */
        int _navdms;
        /**
         * @brief ---
         */
        float _narrow;
        /**
         * @brief ---
         */
        float _aggression;
        /**
         * @brief ---
         */
        int _nsearch;
        /**
         * @brief ---
         */
        float _power;
        /**
         * @brief The threshold for single pulse detection, multiple of standard deviation
         */
        float _sigma_cutoff;
        /**
         * @brief ---
         */
        float _sigma_constant;
        /**
         * @brief ---
         */
        float _max_boxcar_width_in_sec;
        /**
         * @brief ---
         */
        float _wide;
        /**
         * @brief The number of dm ranges
         */
        int _range;
        /**
         * @brief An array containing lowest band of each dm range, specified by the user
         */
        float* _user_dm_low;
        /**
         * @brief An array containing lowest band of each dm range, specified by the user
         */
        float* _user_dm_high;
        /**
         * @brief An array containing lowest band of each dm range, specified by the user
         */
        float* _user_dm_step;
        /**
         * @brief ---
         */
        int* _in_bin;
        /**
         * @brief ---
         */
        int* _out_bin;
        // dedispersion strategy
        // todo: move from float* to vector<float>
        /**
         * @brief Value used to make sure that dms from dm_low to dm_high are used
         */
        int _maxshift;
        /**
         * @brief An array containing the lowest bound of each dm range
         */
        float* _dm_low;
        /**
         * @brief An array containing the highest bound of each dm range
         */
        float* _dm_high;
        /**
         * @brief An array containing the step size of each dm range
         */
        float* _dm_step;
        /**
         * @brief An array containing a constant associated with each channel to perform dedispersion algorithm
         */
        float* _dmshifts;
        /**
         * @brief An array containing the number of dms for each range
         */
        int* _ndms;
        /**
         * @brief The maximum number of dm
         */
        int _max_ndms;
        /**
         * @brief The total number of dm
         */
        int _total_ndms;
        /**
         * @brief The highest dm value
         */
        float _max_dm;
        /**
         * @brief The number of time samples required to search for a dm in each dm range
         */
        int** _t_processed;
        /**
         * @brief The number of bits of the input data
         */
        int _nbits;
        /**
         * @brief The number of IF channels
         */
        int _nifs;
        /**
         * @brief ---
         */
        float _tstart;
        /**
         * @brief Time sample value
         */
        float _tsamp;
        /**
         * @brief The number of time samples
         */
        int _nsamp;
        /**
         * @brief ---
         */
        int _nsamples;
        /**
         * @brief
         */
        int _max_samps;
        /**
         * @brief The number of frequency channels
         */
        int _nchans;
        /**
         * @brief The number of chunks the data are divided in
         */
        unsigned int _num_tchunks;
        /*
         * Frequencies (MHz)
         *
         */
        std::vector<float> _bin_frequencies;
};

} // namespace astroaccelerate


#endif // ASTROACCELERATE_DEDISPERSIONSTRATEGY_H
