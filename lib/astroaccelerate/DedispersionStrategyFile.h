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
     * @details This object fills the Dedispersion Strategy based on an input file.
     * 			This class is intended for users who want to analyse data based on filterbank files as they would do with the C/Cuda master branch
     *
     */
class DedispersionStrategyFile
{
    public:
        /**
        *  @brief Constructor
        *  @details Computes dedispersion strategy based on an input file
        */
        DedispersionStrategyFile(FILE** fp
								,int argc
								,char *argv[]
								,DedispersionStrategy &dedispersion_strategy
								,size_t gpu_memory);

        /**
        *  @brief Destructor
        */
        ~DedispersionStrategyFile();

    private:
        /**
        *  @brief   Get the user input
        *  @details This function read the user input file and stores it. It essentially mirrors the master branch behaviour
        *
        */
        void get_user_input(FILE** fp
        					,int argc
        					,char *argv[]
        					,DedispersionStrategy &dedispersion_strategy);
        /**
         * @brief Read telescope parameters from the header of the input file. It essentially mirrors the master branch behaviour
         */
        void get_file_data(FILE **
        				  ,DedispersionStrategy &DedispersionStrategy);
};

} // namespace astroaccelerate


#endif // ASTROACCELERATE_DEDISPERSIONSTRATEGYFILE_H
