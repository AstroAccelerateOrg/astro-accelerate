#ifndef ASTROACCELERATE_SPS_DEDISPERSIONSTRATEGYFILE_H
#define ASTROACCELERATE_SPS_DEDISPERSIONSTRATEGYFILE_H

#include "../headers/params.h"
#include "../headers/headers_mains.h"
#include "../headers/host_help.h"
#include "DedispersionStrategy.h"

#include <stdio.h>

namespace astroaccelerate {

    /**
     * @brief  Dedispersion Strategy File
     *
     * @details This object files Dedispersion Strategy based on an input file
     *
     */
class DedispersionStrategyFile
{
    public:
        /**
        *  @brief Default constructor
        */
        DedispersionStrategyFile(FILE** fp
								,int argc
								,char *argv[]
								,DedispersionStrategy &dedispersion_strategy);

        /**
        *  @brief Destructor
        */
        ~DedispersionStrategyFile();

    private:
        /**
        *  @brief   Get the user input
        *  @details This function read the user input file and stores it
        *
        */
        void get_user_input(FILE** fp
        					,int argc
        					,char *argv[]
        					,DedispersionStrategy &dedispersion_strategy);
        /**
         * @brief Read telescope parameters from the header of the input file
         */
        void get_file_data(FILE **
        				  ,DedispersionStrategy &DedispersionStrategy);
};

} // namespace astroaccelerate


#endif // ASTROACCELERATE_DEDISPERSIONSTRATEGYFILE_H
