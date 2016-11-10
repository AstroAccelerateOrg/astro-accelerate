#ifndef SKA_ASTROACCELERATE_INPUTDATA_H
#define SKA_ASTROACCELERATE_INPUTDATA_H

#include <stdio.h>

namespace ska {
namespace astroaccelerate {

/**
 * @brief 	Input Data
 * 
 * @details This object carries the input data for both Cpu and Gpu
 * 
 */

class InputData
{
    public:
        /**
        *  @brief Default constructor
        */
        InputData();
        /**
        *  @brief Destructor 
        */
        ~InputData();

    private:
    	// Cpu
    	/**
    	*  @brief 
    	*/
    	unsigned short *_input_buffer;
    	/**
    	*  @brief 
    	*/
    	int 			_nchans;
		/**
    	*  @brief 
    	*/
    	int 			_nsamp;
		/**
    	*  @brief 
    	*/
    	size_t 			_inputsize;
		/**
    	*  @brief 
    	*/
    	int 			**_t_processed;
		/**
    	*  @brief 
    	*/
    	int 			_maxshift;
		// Gpu
		/**
    	*  @brief 
    	*/
    	size_t			_gpu_inputsize;
		/**
    	*  @brief 
    	*/
    	unsigned short 	*_d_input;


};

} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_INPUTDATA_H
