#include "../IOData.h"

namespace ska {
namespace astroaccelerate {
namespace sps{

IOData::IOData()
{
	_input_size 		= 0;
	_input_buffer	= NULL;
	_gpu_input_size 	= 0;
	_d_input 		= NULL;
	_output_size 		= 0;
	_output_buffer	= NULL;
	_gpu_output_size 	= 0;
	_d_output 		= NULL;
}

IOData::~IOData()
{
	// free all the pointers
	free(_input_buffer);
	cudaFree(_d_input);
	// note: it's a float ***, should first free all the _output_buffer[i][j], then
	// the _output_buffer[i], then _output_buffer
	free(_output_buffer);
	cudaFree(_d_output);
}

// setters

void IOData::set_input_size(std::size_t input_size)
{
	_input_size = input_size;
}
void IOData::set_input_buffer(unsigned short* input_buffer)
{
	_input_buffer = input_buffer;
}
void IOData::set_gpu_input_size(std::size_t gpu_input_size)
{
	_gpu_input_size = gpu_input_size;
}
void IOData::set_d_input(unsigned short* d_input)
{
	_d_input = d_input;
}
void IOData::set_output_size(std::size_t output_size)
{
	_output_size = output_size;
}
void IOData::set_output_buffer(float*** output_buffer)
{
	_output_buffer = output_buffer;
}
void IOData::set_gpu_output_size(std::size_t gpu_output_size)
{
	_gpu_output_size = gpu_output_size;
}
void IOData::set_d_output(float* d_output)
{
	_d_output = d_output;
}

// getters
std::size_t IOData::get_input_size() const
{
	return _input_size;
}
unsigned short* IOData::get_input_buffer() const
{
	return _input_buffer;
}
std::size_t IOData::get_gpu_input_size() const
{
	return _gpu_input_size;
}
unsigned short* IOData::get_d_input() const
{
	return _d_input;
}
std::size_t IOData::get_output_size() const
{
	return _output_size;
}
float*** IOData::get_output_buffer() const
{
	return _output_buffer;
}
std::size_t IOData::get_gpu_output_size() const
{
	return _gpu_output_size;
}
float* IOData::get_d_output() const
{
	return _d_output;
}

// methods

void IOData::allocate_memory_cpu_input(DedispersionPlan const &dedispersion_plan)
{
	_input_size = dedispersion_plan.get_nsamp() * (size_t)(dedispersion_plan.get_nchans()) * sizeof(unsigned short);
	_input_buffer = (unsigned short *) malloc(_input_size);
}
void IOData::allocate_memory_cpu_output(DedispersionPlan const &dedispersion_plan)
{
	_output_size = 0;
	_output_buffer = (float ***) malloc(dedispersion_plan.get_range() * sizeof(float **));
	for (int i = 0; i < dedispersion_plan.get_range(); ++i)
	{
		int total_samps = 0;
		for (int k = 0; k < dedispersion_plan.get_num_tchunks(); ++k)
			total_samps += dedispersion_plan.get_t_processed()[i][k];
		_output_buffer[i] = (float **) malloc(dedispersion_plan.get_ndms()[i] * sizeof(float *));
		for (int j = 0; j < dedispersion_plan.get_ndms()[i]; ++j)
		{
			_output_buffer[i][j] = (float *) malloc(( total_samps ) * sizeof(float));
		}
		_output_size += ( total_samps ) * dedispersion_plan.get_ndms()[i] * sizeof(float);
	}
}
void IOData::allocate_memory_gpu(DedispersionPlan const &dedispersion_plan)
{

	int time_samps = dedispersion_plan.get_t_processed()[0][0] + dedispersion_plan.get_maxshift();
	_gpu_input_size = time_samps * (size_t) dedispersion_plan.get_nchans() * sizeof(unsigned short);
	( cudaMalloc((void **) _d_input, _gpu_input_size) );

	if (dedispersion_plan.get_nchans() < dedispersion_plan.get_max_ndms())
	{
		_gpu_output_size = time_samps * dedispersion_plan.get_max_ndms() * sizeof(float);
	}
	else
	{
		_gpu_output_size = time_samps * dedispersion_plan.get_nchans() * sizeof(float);
	}
	( cudaMalloc((void **) _d_output, _gpu_output_size) );

	( cudaMemset(_d_output, 0, _gpu_output_size) );
}

void IOData::get_recorded_data(FILE **fp, DedispersionPlan const &dedispersion_plan)
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
