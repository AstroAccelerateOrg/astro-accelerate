#include "../Sps.h"

#include "../../AstroAccelerate/device_init.h"
#include "../../AstroAccelerate/device_load_data.h"
#include "../../AstroAccelerate/device_zero_dm.h"
#include "../../AstroAccelerate/device_corner_turn.h"
#include "../../AstroAccelerate/device_dedisperse.h"
#include "../../AstroAccelerate/device_bin.h"
#include "../../AstroAccelerate/device_save_data.h"
#include "../../AstroAccelerate/host_write_file.h"
#include "../../AstroAccelerate/host_analysis.h"

namespace ska {
namespace astroaccelerate {
namespace sps {

template<typename SpsParameterType>
Sps<SpsParameterType>::Sps(IOData &io_data,
						   DedispersionPlan  &dedispersion_plan,
						   UserInput &user_input)
						   : _num_tchunks(dedispersion_plan.get_num_tchunks()),
						     _range(dedispersion_plan.get_range()),
						     _nchans(dedispersion_plan.get_nchans()),
						     _maxshift(dedispersion_plan.get_maxshift()),
						     _tsamp(dedispersion_plan.get_tsamp()),
						     _max_ndms(dedispersion_plan.get_max_ndms()),
						     _sigma_cutoff(user_input.get_sigma_cutoff()),
						     _gpu_outputsize(io_data.get_gpu_output_size()),
						     _ndms(dedispersion_plan.get_ndms()),
						     _dmshifts(dedispersion_plan.get_dmshifts()),
						     _t_processed(dedispersion_plan.get_t_processed()),
						     _dm_low(dedispersion_plan.get_dm_low()),
						     _dm_high(dedispersion_plan.get_dm_high()),
						     _dm_step(dedispersion_plan.get_dm_step()),
						     _in_bin(user_input.get_in_bin()),
						     _out_bin(user_input.get_out_bin()),
						     _input_buffer(io_data.get_input_buffer()),
						     _output_buffer(io_data.get_output_buffer()),
						     _d_input(io_data.get_d_input()),
						     _d_output(io_data.get_d_output())
{
}

template<typename SpsParameterType>
Sps<SpsParameterType>::~Sps()
{
}

template<typename SpsParameterType>
void Sps<SpsParameterType>::operator()( unsigned device_id,
										IOData &io_data,
										DedispersionPlan &dedispersion_plan,
                                        UserInput &user_input)
{
	//
		long int inc = 0;
		float tstart_local = 0.0f;

		// Initialise the GPU.
		size_t gpu_memory = 0;
		cudaSetDevice(device_id);
		size_t mem_free, total;
		cudaMemGetInfo(&mem_free, &total);
		gpu_memory = ( mem_free/4 );

		// Call the strategy method of class dedispersion (could be done outside ?)
		// -> can't be done outside b/c need gpu_mem and gpu is initialised inside the library
		// Call the strategy method
		dedispersion_plan.make_strategy(user_input.get_user_dm_low(),
										user_input.get_user_dm_high(),
										user_input.get_user_dm_step(),
										user_input.get_in_bin(),
										gpu_memory
										);
		// allocate memory cpu output
		io_data.allocate_memory_cpu_output(dedispersion_plan);
		// allocate memory gpu
		io_data.allocate_memory_gpu(dedispersion_plan);

		//printf("De-dispersing...");
		int t, dm_range;
		float tsamp_original = _tsamp;
		int maxshift_original = _maxshift;

		//
		float *out_tmp = NULL;
		out_tmp = (float *) malloc(( _t_processed[0][0] + _maxshift ) * _max_ndms * sizeof(float));
		memset(out_tmp, 0.0f, _t_processed[0][0] + _maxshift * _max_ndms * sizeof(float));

		for (t = 0; t < _num_tchunks; ++t)
		{

			load_data(-1, _in_bin, _d_input, &_input_buffer[(long int) ( inc * _nchans )],
								_t_processed[0][t], _maxshift, _nchans, _dmshifts);

			if (user_input.get_enable_zero_dm())
				zero_dm(_d_input, _nchans, _t_processed[0][t]+_maxshift);

			corner_turn(_d_input, _d_output, _nchans, _t_processed[0][t] + _maxshift);
			int oldBin = 1;

			for (dm_range = 0; dm_range < _range; ++dm_range)
			{

				_maxshift = maxshift_original / _in_bin[dm_range];

				cudaDeviceSynchronize();
				load_data(dm_range, _in_bin, _d_input, &_input_buffer[(long int) ( inc * _nchans )], _t_processed[dm_range][t], _maxshift, _nchans, _dmshifts);

				if (_in_bin[dm_range] > oldBin)
				{
					bin_gpu(_d_input, _d_output, _nchans, _t_processed[dm_range - 1][t] + _maxshift * _in_bin[dm_range]);
					( _tsamp ) = ( _tsamp ) * 2.0f;
				}

				dedisperse(dm_range, _t_processed[dm_range][t], _in_bin, _dmshifts, _d_input, _d_output, _nchans,
				           ( _t_processed[dm_range][t] + _maxshift ), _maxshift, &_tsamp, _dm_low, _dm_high, _dm_step, _ndms);

				_gpu_outputsize = _ndms[dm_range] * ( _t_processed[dm_range][t] ) * sizeof(float);

				save_data(_d_output, out_tmp, _gpu_outputsize);

				//#pragma omp parallel for
				for (int k = 0; k < _ndms[dm_range]; ++k)
				{
					memcpy(&_output_buffer[dm_range][k][inc / _in_bin[dm_range]], &out_tmp[k * _t_processed[dm_range][t]],
								sizeof(float) * _t_processed[dm_range][t]);
				}

				if (user_input.get_output_dmt() == 1)
					write_output(dm_range, _t_processed[dm_range][t], _ndms[dm_range], gpu_memory, out_tmp, _gpu_outputsize, _dm_low, _dm_high);
				if (user_input.get_enable_analysis() == 1)
					analysis(dm_range, tstart_local, _t_processed[dm_range][t], ( _t_processed[dm_range][t] + _maxshift ), _nchans, _maxshift, _max_ndms, _ndms, _out_bin, _sigma_cutoff, _d_output, _dm_low, _dm_high, _dm_step, _tsamp);
				oldBin = _in_bin[dm_range];
			}

			memset(out_tmp, 0.0f, _t_processed[0][0] + _maxshift * _max_ndms * sizeof(float));

			inc = inc + _t_processed[0][t];
			printf("\nINC:\t%ld", inc);
			tstart_local = ( tsamp_original * inc );
			_tsamp = tsamp_original;
			_maxshift = maxshift_original;
		}

		free(out_tmp);
}

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
