#include "../InputData.h"

namespace ska {
namespace astroaccelerate {
namespace sps{

InputData::InputData()
{
	_input_size 		= 0;
	_input_buffer	= NULL;
}

InputData::~InputData()
{
	// free all the pointers
	free(_input_buffer);
}

// setters

void InputData::set_input_size(std::size_t input_size)
{
	_input_size = input_size;
}
void InputData::set_input_buffer(unsigned short* input_buffer)
{
	_input_buffer = input_buffer;
}

// getters
std::size_t InputData::get_input_size() const
{
	return _input_size;
}
unsigned short* InputData::get_input_buffer() const
{
	return _input_buffer;
}

// methods

void InputData::allocate_memory_cpu_input(DedispersionPlan const &dedispersion_plan)
{
	_input_size = dedispersion_plan.get_nsamp() * (size_t)(dedispersion_plan.get_nchans()) * sizeof(unsigned short);
	_input_buffer = (unsigned short *) malloc(_input_size);
}

void InputData::get_recorded_data(FILE **fp, DedispersionPlan const &dedispersion_plan)
{

	int c;
	unsigned long int total_data;
	//{{{ Load in the raw data from the input file and transpose
	if (dedispersion_plan.get_nbits() == 8)
	{
		// Allocate a tempory buffer to store a line of frequency data
		unsigned char *temp_buffer = (unsigned char *) malloc(dedispersion_plan.get_nchans() * sizeof(unsigned char));
		// Read in the data, transpose it and store it in the input buffer
		total_data = 0;
		while (!feof(*fp))
		{
			if (fread(temp_buffer, sizeof(unsigned char), dedispersion_plan.get_nchans(), *fp) != dedispersion_plan.get_nchans())
				break;
			for (c = 0; c < dedispersion_plan.get_nchans(); c++)
			{
				_input_buffer[c + total_data * ( dedispersion_plan.get_nchans() )] = (unsigned short) temp_buffer[c];
			}
			total_data++;

		}
		free(temp_buffer);
	}
	else
	{
		printf("\n\n========================= ERROR =========================\n");
		printf(" This is a SKA prototype code and only runs with 8 bit data\n");
		printf("\n=========================================================\n");
	}
	//}}}
}

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
