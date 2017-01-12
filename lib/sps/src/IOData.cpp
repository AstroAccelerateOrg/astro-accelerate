#include "../IOData.h"

namespace ska {
namespace astroaccelerate {
namespace sps{

IOData::IOData()
{
	_input_size 		= 0;
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
unsigned short* IOData::get_input_buffer()
{
	return _input_buffer.data();
}
std::size_t IOData::get_gpu_input_size() const
{
	return _gpu_input_size;
}
unsigned short* IOData::get_d_input()
{
	return _d_input;
}
std::size_t IOData::get_output_size() const
{
	return _output_size;
}
float*** IOData::get_output_buffer()
{
	return _output_buffer;
}
std::size_t IOData::get_gpu_output_size() const
{
	return _gpu_output_size;
}
float* IOData::get_d_output()
{
	return _d_output;
}

// methods

void IOData::allocate_memory_cpu_input(DedispersionPlan &dedispersion_plan)
{
	_input_size = dedispersion_plan.get_nsamp() * (size_t)(dedispersion_plan.get_nchans()) * sizeof(unsigned short);
	_input_buffer.resize(_input_size);
}

void IOData::allocate_memory_cpu_output(DedispersionPlan &dedispersion_plan)
{
	int range 				= dedispersion_plan.get_range();
	int num_tchunks 	= dedispersion_plan.get_num_tchunks();
	int** t_processed = dedispersion_plan.get_t_processed();
	int* ndms 				= dedispersion_plan.get_ndms();

	_output_size = 0;
	_output_buffer = (float ***) malloc(range * sizeof(float **));
	for (int i = 0; i < range; ++i)
	{
		int total_samps = 0;
		for (int k = 0; k < num_tchunks; ++k)
			total_samps += t_processed[i][k];
		//printf("\nTOTSAMPS:\t%d %d", total_samps, i);
		_output_buffer[i] = (float **) malloc(ndms[i] * sizeof(float *));
		//if((*output_buffer)[i]) printf("\n FAILED! Could not allocate %zu bytes", ndms[i]*sizeof(float *));
		for (int j = 0; j < ndms[i]; ++j)
		{
			_output_buffer[i][j] = (float *) malloc(( total_samps ) * sizeof(float));
			//if((*output_buffer)[i][j]) printf("\n FAILED! Could not allocate %zu bytes", ndms[i]*sizeof(float *));
//			memset((*output_buffer)[i][j],0.0f,(total_samps)*sizeof(float));
		}
		_output_size += ( total_samps ) * ndms[i] * sizeof(float);
		printf("\noutput size: %llu", (unsigned long long) sizeof( _output_buffer ) / 1024 / 1024 / 1024);
	}
}

void IOData::allocate_memory_gpu(DedispersionPlan &dedispersion_plan)
{
	int maxshift = dedispersion_plan.get_maxshift();
	int** t_processed = dedispersion_plan.get_t_processed();
	int nchans = dedispersion_plan.get_nchans();
	int max_ndms = dedispersion_plan.get_max_ndms();

	int time_samps = t_processed[0][0] + maxshift;
	_gpu_input_size = time_samps * (size_t) nchans * sizeof(unsigned short);
	( cudaMalloc((void **) _d_input, _gpu_input_size) );

	if (nchans < max_ndms)
	{
		_gpu_output_size = time_samps * max_ndms * sizeof(float);
	}
	else
	{
		_gpu_output_size = time_samps * nchans * sizeof(float);
	}
	( cudaMalloc((void **) _d_output, _gpu_output_size) );

	//end_t=omp_get_wtime();
	//time = (float)(end_t-start_t);
	//printf("\nGPU Malloc in: %f ", time);

	( cudaMemset(_d_output, 0, _gpu_output_size) );
}

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
