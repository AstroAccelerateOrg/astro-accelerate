#include "AstroAccelerate/OutputData.h"

namespace ska {
namespace astroaccelerate {

OutputData::OutputData()
{
	// Cpu
    _output_buffer	= NULL;
    _outputsize 	= 0;
    _range  		= 0;
    _t_processed = NULL;	// = NULL looks better ?
    // Gpu
    _d_output 		= NULL;
    _max_ndms 		= 0;
    _maxshift 		= 0;
    _nchans			= 0;
    _gpu_outputsize	= 0;
}

OutputData::~OutputData()
{
}

} // namespace astroaccelerate
} // namespace ska
