#ifndef SKA_ASTROACCELERATE_SPS_OUTPUTDATA_H
#define SKA_ASTROACCELERATE_SPS_OUTPUTDATA_H

#include <stdio.h>

namespace ska {
namespace astroaccelerate {

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
    	/**
    	*  @brief 
    	*/
    	int		_range;
    	/**
    	*  @brief 
    	*/
    	int		**_t_processed;
    	// Gpu
    	/**
    	*  @brief 
    	*/
    	float 	*_d_output;
    	/**
    	*  @brief 
    	*/
    	int		_max_ndms;
    	/**
    	*  @brief 
    	*/
    	int 	_maxshift;
    	/**
    	*  @brief 
    	*/
    	int 	_nchans;
    	/**
    	*  @brief 
    	*/
    	size_t 	_gpu_outputsize;
};

} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_OUTPUTDATA_H
