#include "../InputData.h"


namespace ska {
namespace astroaccelerate {
namespace sps{


InputData::InputData()
{
	// Cpu
	_input_size 		= 0;
	_input_buffer	= NULL;
	// Gpu
	_gpu_input_size 	= 0;
	_d_input 		= NULL;
}

InputData::~InputData()
{
}

void InputData::set_input_size(size_t input_size)
{
	_input_size = input_size;
}

void InputData::set_input_buffer(unsigned short* input_buffer)
{
	_input_buffer = input_buffer;
}

void InputData::set_gpu_input_size(size_t gpu_input_size)
{
	_gpu_input_size = gpu_input_size;
}

void InputData::set_input_buffer(unsigned short* input_buffer)
{
	_d_input = d_input;
}

size_t InputData::get_input_size() const
{
	return _input_size;
}
unsigned short* InputData::get_input_buffer() const
{
	return _input_buffer;
}
size_t InputData::get_gpu_input_size() const
{
	return _gpu_input_size;
}
unsigned short* InputData::get_d_input() const
{
	return _d_input;
}

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
