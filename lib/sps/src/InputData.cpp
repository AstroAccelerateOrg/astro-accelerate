#include "AstroAccelerate/InputData.h"

namespace ska {
namespace astroaccelerate {


InputData::InputData()
{



	// Cpu
	_inputsize 		= 0;
	_input_buffer	= NULL;
	// Gpu
	_gpu_inputsize 	= 0;
	_d_input 		= NULL;
}

InputData::~InputData()
{
}

} // namespace astroaccelerate
} // namespace ska
