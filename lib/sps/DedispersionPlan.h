#ifndef SKA_ASTROACCELERATE_SPS_DEDISPERSIONPLAN_H
#define SKA_ASTROACCELERATE_SPS_DEDISPERSIONPLAN_H

#include <stdio.h>
#include "UserInput.h"
#include <vector>

namespace ska {
namespace astroaccelerate {
namespace sps {

	/**
	 * @brief 	Dedispersion Plan
	 *
	 * @details	This object carries the dedispersion plan
	 *
	 */

class DedispersionPlan
{
    public:
        /**
        *  @brief Default constructor
        */
        DedispersionPlan();

        /**
        *  @brief Destructor
        */
        ~DedispersionPlan();

        void add_range(float begin, float end, float step);
        // add 1 to in and out_bin

        /**
	       *  @brief Setters
	      */
        void set_in_bin(int *);
        void set_out_bin(int *);
        void set_maxshift(int);
        void set_max_ndms(int);
        void set_total_ndms(int);
        void set_range(int);
        void set_nbits(int);
        void set_nifs(int);
        void set_tstart(float);
        void set_tsamp(float);
        void set_nsamp(int);
        void set_nsamples(int);
        void set_max_samps(int);
        void set_nchans(int);
        void set_fch1(float);
        void set_foff(float);
        void set_num_tchunks(unsigned int);
        void set_power(float);

        /**
        *  @brief Getters
        */
        int*   get_in_bin();
        int*   get_out_bin();
        int    get_maxshift() const;
        float* get_dm_low();
        float* get_dm_high();
        float* get_dm_step();
        float* get_dmshifts();
        int*   get_ndms();
        int 		 get_max_ndms() const;
        int 		 get_total_ndms() const;
        float 		 get_max_dm() const;
        int  		 get_range() const;
        int**  		 get_t_processed();
        int 		 get_nbits() const;
        int 		 get_nifs() const;
        float 		 get_tstart() const;
        float  		 get_tsamp() const;
        int  		 get_nsamp() const;
        int			 get_nsamples() const;
        int 		 get_max_samps() const;
        int  		 get_nchans() const;
        float  		 get_fch1() const;
        float  		 get_foff() const;
      	unsigned int get_num_tchunks() const;
      	float	     get_power() const;

        /**
         * @brief Return the minimum number of samples required for the
         *        algorithm to operate
         */
        int minimum_number_of_samples() const;

    private:
        void make_strategy(float*, float*, float*, size_t);

    private:
      int* 		_in_bin;
      int* 		_out_bin;
      int 		_maxshift;
      std::vector<float> 	_dm_low;
      std::vector<float> 	_dm_high;
      std::vector<float> 	_dm_step;
      std::vector<float> 	_dmshifts;
      std::vector<int> 		_ndms;
    	int 		_max_ndms;
    	int 		_total_ndms;
    	float 	_max_dm;
    	int 		_range;
    	int** 	_t_processed;
    	int 		_nbits;
    	int 		_nifs;
    	float 	_tstart;
    	float 	_tsamp;
    	int 		_nsamp;
    	int 		_nsamples;
    	int 		_max_samps;
    	int 		_nchans;
    	float 	_fch1;
    	float 	_foff;
    	unsigned int 		_num_tchunks;
    	float 	_power;
    	// add num_tchunks
};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_DEDISPERSIONPLAN_H
