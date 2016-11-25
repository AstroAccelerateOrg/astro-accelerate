#ifndef SKA_ASTROACCELERATE_USERINPUT_H
#define SKA_ASTROACCELERATE_USERINPUT_H

#include <stdio.h>

namespace ska {
namespace astroaccelerate {

/**
 * @brief
 * 
 * @details
 * 
 */

class UserInput
{
    public:
				/**
				*  @brief Constructor
				*/
				UserInput();
        /**
        *  @brief Destructor
        */
        ~UserInput();
        /**
        *  @brief Getters
        */
        int 		get_multi_file() const;
        int 		get_enable_debug() const;
        int 		get_enable_analysis() const;
        int 		get_enable_periodicity() const;
        int 		get_enable_acceleration() const;
        int 		get_output_dm() const;
        int 		get_enable_zero_dm() const;
        int 		get_nboots() const ;
        int 		get_ntrial_bons() const;
        int 		get_navdms() const;
        int 		get_narrow() const;
        int 		get_agression() const;
        int			get_nsearch() const;
        int* 		get_inBin() const;
        int* 		get_outBin() const;
        float 	get_power() const;
        float 	get_sigma_cutoff() const;
        int 		get_range() const;
        float*	get_user_dm_low() const;
        float*	get_user_dm_high() const;
        float* 	get_user_dm_step() const;

        /**
        *  @brief Get the user input
        *  				This is not a class Getter though
        *  				This method replaces all the Setters
        *  				Name may be confusing
        */
        void 	get_user_input(FILE** fp, int argc, char *argv[]);

    private:
        int 		_multi_file;
        int 		_enable_debug;
        int 		_enable_analysis;
        int 		_enable_periodicity;
        int			_enable_acceleration;
        int 		_output_dm;
        int 		_enable_zero_dm;
        int 		_nboots;
        int 		_ntrial_bins;
        int 		_navdms;
        int 		_narrow;
        int 		_agression;
        int			_nsearch;
        int*		_inBin;
        int* 		_outBin;
        float 	_power;
        float 	_sigma_cutoff;
        int 		_range;
        float* 	_user_dm_low;
        float* 	_user_dm_high;
        float*	_user_dm_step;
};

} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_USERINPUT_H
