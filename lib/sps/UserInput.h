#ifndef SKA_ASTROACCELERATE_SPS_USERINPUT_H
#define SKA_ASTROACCELERATE_SPS_USERINPUT_H

#include <stdio.h>

namespace ska {
namespace astroaccelerate {
namespace sps{

/**
 * @brief   User Input
 * 
 * @details This object read the user input file and carries its informations
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
        int     get_multi_file() const;
        int     get_enable_debug() const;
        int     get_enable_analysis() const;
        int     get_enable_periodicity() const;
        int     get_enable_acceleration() const;
        int     get_output_dmt() const;
        int     get_enable_zero_dm() const;
        int     get_nboots() const ;
        int     get_ntrial_bins() const;
        int     get_navdms() const;
        float   get_narrow() const;
        float   get_aggression() const;
        int     get_nsearch() const;
        float   get_power() const;
        float   get_sigma_cutoff() const;
        float   get_wide() const;
        int     get_range() const;
        float*  get_user_dm_low() const;
        float*  get_user_dm_high() const;
        float*  get_user_dm_step() const;
        int*    get_in_bin() const;
        int*    get_out_bin() const;

        /**
        *  @brief   Get the user input
        *  @details This function read the user input file and stores it
        *
        */
        void    get_user_input(FILE** fp, int argc, char *argv[]);

    private:
        int     _multi_file;
        int     _enable_debug;
        int     _enable_analysis;
        int     _enable_periodicity;
        int     _enable_acceleration;
        int     _output_dmt;
        int     _enable_zero_dm;
        int     _nboots;
        int     _ntrial_bins;
        int     _navdms;
        float   _narrow;
        float   _aggression;
        int     _nsearch;
        float   _power;
        float   _sigma_cutoff;
        float   _wide;
        int     _range;
        float*  _user_dm_low;
        float*  _user_dm_high;
        float*  _user_dm_step;
        int*    _in_bin;
        int*    _out_bin;
};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_SPS_USERINPUT_H
