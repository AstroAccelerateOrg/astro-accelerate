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

        /**
         * @brief Return the minimum number of samples required for the
         *        algorithm to operate
         */
        int minimum_number_of_samples() const;

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
};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_SPS_DEDISPERSIONPLAN_H
