#include "AstroAccelerate/InputData.h"

namespace ska {
namespace astroaccelerate {


InputData::InputData()
{



	// Cpu
	_input_buffer	= NULL;
    _nchans			= 0;
	_nsamp			= 0;
	_inputsize 		= 0;
	_t_processed = NULL;
	_maxshift		= 0;
	// Gpu
	_inputsize 		= 0;
	_gpu_inputsize 	= 0;
	_d_input 		= NULL;
}

InputData::~InputData()
{
}

} // namespace astroaccelerate
} // namespace ska
