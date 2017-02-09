#ifndef SKA_ASTROACCELERATE_SPS_INPUTDATA_H
#define SKA_ASTROACCELERATE_SPS_INPUTDATA_H

//
#include <driver_functions.h>
#include <vector_types.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cstdlib>

//
#include "DedispersionPlan.h"

namespace ska {
namespace astroaccelerate {
namespace sps{

/**
 * @brief 	Input/Output Data
 *
 * @details This object carries the input/output data of SPS for both Cpu and Gpu
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
                void set_input_size(std::size_t);
                void set_input_buffer(unsigned short*);
                /**
                *  @brief Getters
                */
                std::size_t get_input_size() const;
                unsigned short* get_input_buffer() const;
                /**
                 *
                 */
                void allocate_memory_cpu_input(DedispersionPlan const &);
                void get_recorded_data(FILE **, DedispersionPlan  const &);

    private:
                /**
                *  @brief   Host/Cpu input members
                *  @details _input_buffer carries _input_size of data on the host
                *  @see allocate_memory_cpu_input()
                */
                size_t          _input_size;
                unsigned short *_input_buffer;
};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_INPUTDATA_H
