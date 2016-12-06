#ifndef SKA_ASTROACCELERATE_SPS_INPUTDATA_H
#define SKA_ASTROACCELERATE_SPS_INPUTDATA_H

#include <stdio.h>

namespace ska {
namespace astroaccelerate {
namespace sps{


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
        /**
        *  @brief Setters
        */
        void set_input_size(size_t);
        void set_input_buffer(unsigned short*);
        void set_gpu_input_size(size_t);
        void set_d_input(unsigned short*);
        /**
        *  @brief Getters
        */
        size_t get_input_size() const;
        unsigned short* get_input_buffer() const;
        size_t get_gpu_input_size() const;
        unsigned short* get_d_input() const;


    private:
    	// Cpu
    	size_t 			_input_size;
    	unsigned short *_input_buffer;
		// Gpu
    	size_t			_gpu_input_size;
    	unsigned short 	*_d_input;

};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_INPUTDATA_H
