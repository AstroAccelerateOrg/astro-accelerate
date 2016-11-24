#ifndef SKA_ASTROACCELERATE_SPS_DEDISPERSIONPLAN_H
#define SKA_ASTROACCELERATE_SPS_DEDISPERSIONPLAN_H

#include <stdio.h>

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

        /**
	       *  @brief Setters
	      */
        void 		set_in_bin(int *);
        void 		set_out_bin(int *);
        void 		set_maxshift(int);
	      void    set_user_dm_low(float *);
	      void    set_user_dm_high(float *);
	      void    set_user_dm_step(float *);
        void    set_dm_low(float *);
        void    set_dm_high(float *);
        void    set_dm_step(float *);
        void    set_dmshifts(float *);
        void    set_range(int);
        void    set_ndms(int *);
        void    set_t_processed(int **);
        void    set_tsamp(float);
        void    set_max_ndms(int);
        void    set_nsamp(int);
        void    set_nchans(int);
        void    set_fch1(float);
        void    set_foff(float);

        /**
        *  @brief Getters
        */
        int*		get_in_bin() const;
        int* 		get_out_bin() const;
        int 		get_maxshift() const;
        float*  get_user_dm_low() const;
        float*  get_user_dm_high() const;
        float*  get_user_dm_step() const;
        float*  get_dm_low() const;
        float*  get_dm_high() const;
        float*  get_dm_step() const;
        float*  get_dmshifts() const;

        int  		get_range() const;
        int*  	get_ndms() const;
        int**  	get_t_processed() const;
        float  	get_tsamp() const;
        int  		get_max_ndms() const;
        int  		get_nsamp() const;
        int  		get_nchans() const;
        float  	get_fch1() const;
        float  	get_foff() const;

        /**
         * @brief Return the minimum number of samples required for the
         *        algorithm to operate
         */
        int minimum_number_of_samples() const;

    private:
        void make_strategy(float* user_dm_low, float* user_dm_high, float* user_dm_steps);

    private:
      /**
      *  @brief
      */
      int *_in_bin;
      /**
      *  @brief
      */
      int *_out_bin;
      /**
      *  @brief
      */
      int _maxshift;
      /**
    	*  @brief
    	*/
    	float *_user_dm_low;
			/**
    	*  @brief
    	*/
    	float *_user_dm_high;
			/**
    	*  @brief
    	*/
    	float *_user_dm_step;
			/**
    	*  @brief
    	*/
    	float *_dm_low;
			/**
    	*  @brief
    	*/
    	float *_dm_high;
			/**
    	*  @brief
    	*/
    	float *_dm_step;
			/**
    	*  @brief
    	*/
    	float *_dmshifts;
    	/**
    	*  @brief
    	*/
    	int _range;
    	int *_ndms;
    	int **_t_processed;
    	float _tsamp;
    	int _max_ndms;
    	int _nsamp;
    	int _nchans;
    	float _fch1;
    	float _foff;

};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_SPS_DEDISPERSIONPLAN_H
