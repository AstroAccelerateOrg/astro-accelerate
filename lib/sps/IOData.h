#ifndef SKA_ASTROACCELERATE_SPS_IODATA_H
#define SKA_ASTROACCELERATE_SPS_IODATA_H

/* This function takes a pointer to the file pointer so that it can update the position of the file pointer
 */
#include <vector_types.h>
#include <driver_functions.h>
#include <cuda_runtime.h>

// CUDA utilities and system includes
#include <vector_types.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <omp.h>
#include <cuda.h>

#include <cstdlib>
#include <vector>

#include "DedispersionPlan.h"

namespace ska {
namespace astroaccelerate {
namespace sps{

/**
 * @brief 	Input/Output Data
 * 
 * @details This object carries the input/output data for both Cpu and Gpu
 * 
 */

class IOData
{
    public:
        /**
        *  @brief Default constructor
        */
        IOData();
        /**
        *  @brief Destructor 
        */
        ~IOData();
        /**
        *  @brief Setters
        */
        void set_input_size(std::size_t);
        void set_gpu_input_size(std::size_t);
        void set_d_input(unsigned short*);
        void set_output_size(std::size_t);
        void set_gpu_output_size(std::size_t);
        void set_d_output(float*);
        /**
        *  @brief Getters
        */
        std::size_t get_input_size() const;
        unsigned short* get_input_buffer();
        std::size_t get_gpu_input_size() const;
        unsigned short* get_d_input();
        std::size_t get_output_size() const;
        float*** get_output_buffer();
        std::size_t get_gpu_output_size() const;
        float* get_d_output();

        void allocate_memory_cpu_input(DedispersionPlan &);
        void allocate_memory_cpu_output(DedispersionPlan &);
        void allocate_memory_gpu(DedispersionPlan &);

    private:
      // INPUT
    	// Cpu
    	size_t 			_input_size;
    	std::vector<unsigned short> _input_buffer;
    	// Gpu
    	size_t			_gpu_input_size;
    	unsigned short 	*_d_input;
    	// OUTPUT
    	// Cpu
    	size_t 	_output_size;
    	float 	***_output_buffer;
    	// Gpu
    	size_t 	_gpu_output_size;
    	float 	*_d_output;

};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_IODATA_H
