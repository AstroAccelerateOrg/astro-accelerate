#ifndef SKA_ASTROACCELERATE_SPS_OUTPUTDATA_H
#define SKA_ASTROACCELERATE_SPS_OUTPUTDATA_H

#include <stdio.h>

namespace ska {
namespace astroaccelerate {
namespace sps{
/**
 * @brief 	Output Data
 * 
 * @details This object carries the output data for both Cpu and Gpu
 * 
 */

class OutputData
{
    public:
        /**
        *  @brief Default constructors
        */
        OutputData();
        /**
        *  @brief Destructor 
        */
        ~OutputData();

    private:
    	// Cpu
    	/**
    	*  @brief 
    	*/
    	float 	***_output_buffer;
    	/**
    	*  @brief 
    	*/
    	size_t 	_outputsize;

    	// Gpu
    	/**
    	*  @brief 
    	*/
    	float 	*_d_output;
    	/**
    	*  @brief 
    	*/
    	size_t 	_gpu_outputsize;
};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_OUTPUTDATA_H
