#ifndef SKA_ASTROACCELERATE_DEDISPERSIONPLAN_H
#define SKA_ASTROACCELERATE_DEDISPERSIONPLAN_H

#include <stdio.h>

namespace ska {
namespace astroaccelerate {

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
        float*  get_user_dm_low();
        float*  get_user_dm_high();
        float*  get_user_dm_step();
        float*  get_dm_low();
        float*  get_dm_high();
        float*  get_dm_step();
        float*  get_dmshifts();
    private:
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

} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_DEDISPERSIONPLAN_H
