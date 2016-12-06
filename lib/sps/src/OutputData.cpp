#include "../OutputData.h"

namespace ska {
namespace astroaccelerate {
namespace sps{


OutputData::OutputData()
{
	// Cpu
    _output_buffer	= NULL;
    _outputsize 	= 0;
    // Gpu
    _d_output 		= NULL;
    _gpu_outputsize	= 0;
}

OutputData::~OutputData()
{
}

}
} // namespace astroaccelerate
} // namespace ska
