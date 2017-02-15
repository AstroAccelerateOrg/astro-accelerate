#include "../Sps.h"

#include "../../AstroAccelerate/headers_mains.h"
#include "../../AstroAccelerate/device_bin.h"
#include "../../AstroAccelerate/device_init.h"
#include "../../AstroAccelerate/device_dedisperse.h"
#include "../../AstroAccelerate/device_dedispersion_kernel.h"
#include "../../AstroAccelerate/device_zero_dm.h"

#include "../../AstroAccelerate/device_SPS_inplace_kernel.h" //Added by KA
#include "../../AstroAccelerate/device_SPS_inplace.h" //Added by KA
#include "../../AstroAccelerate/device_MSD_grid.h" //Added by KA
#include "../../AstroAccelerate/device_MSD_plane.h" //Added by KA
#include "../../AstroAccelerate/device_MSD_limited.h" //Added by KA
#include "../../AstroAccelerate/device_SNR_limited.h" //Added by KA
#include "../../AstroAccelerate/device_threshold.h" //Added by KA
#include "../../AstroAccelerate/device_single_FIR.h" //Added by KA

#include "../../AstroAccelerate/device_load_data.h"
#include "../../AstroAccelerate/device_corner_turn.h"
#include "../../AstroAccelerate/device_save_data.h"
#include "../../AstroAccelerate/host_acceleration.h"
#include "../../AstroAccelerate/host_allocate_memory.h"
#include "../../AstroAccelerate/host_analysis.h"
#include "../../AstroAccelerate/host_periods.h"
#include "../../AstroAccelerate/host_debug.h"
#include "../../AstroAccelerate/host_get_file_data.h"
#include "../../AstroAccelerate/host_get_recorded_data.h"
#include "../../AstroAccelerate/host_get_user_input.h"
#include "../../AstroAccelerate/host_help.h"
#include "../../AstroAccelerate/host_rfi.h"
#include "../../AstroAccelerate/host_stratagy.h"
#include "../../AstroAccelerate/host_write_file.h"

#include "../../AstroAccelerate/params.h"

#include "../../AstroAccelerate/host_debug.h"


namespace ska {
namespace astroaccelerate {
namespace sps {

template<typename SpsParameterType>
Sps<SpsParameterType>::Sps(InputData &input_data,
		 DedispersionPlan  &dedispersion_plan,
		 UserInput &user_input)
		 : _num_tchunks(dedispersion_plan.get_num_tchunks()),
		 _range(dedispersion_plan.get_range()),
		 _nchans(dedispersion_plan.get_nchans()),
		 _maxshift(dedispersion_plan.get_maxshift()),
		 _tsamp(dedispersion_plan.get_tsamp()),
		 _max_ndms(dedispersion_plan.get_max_ndms()),
		 _sigma_cutoff(user_input.get_sigma_cutoff()),
		 _ndms(dedispersion_plan.get_ndms()),
		 _dmshifts(dedispersion_plan.get_dmshifts()),
		 _t_processed(dedispersion_plan.get_t_processed()),
		 _dm_low(dedispersion_plan.get_dm_low()),
		 _dm_high(dedispersion_plan.get_dm_high()),
		 _dm_step(dedispersion_plan.get_dm_step()),
		 _in_bin(user_input.get_in_bin()),
		 _out_bin(user_input.get_out_bin()),
		 _input_buffer(input_data.get_input_buffer())
{
	_gpu_input_size 	= 0;
	_d_input 			= NULL;
	_gpu_output_size 	= 0;
	_d_output 			= NULL;
}

template<typename SpsParameterType>
Sps<SpsParameterType>::~Sps()
{
	cudaFree(_d_input);
	cudaFree(_d_output);
}


template<typename SpsParameterType>
void Sps<SpsParameterType>::allocate_memory_cpu_output(DedispersionPlan const &dedispersion_plan,
													   float ****output_buffer,
													   size_t *output_size)
{
	*output_buffer = (float ***) malloc(dedispersion_plan.get_range() * sizeof(float **));
	for (int i = 0; i < dedispersion_plan.get_range(); ++i)
	{
		int total_samps = 0;
		for (int k = 0; k < dedispersion_plan.get_num_tchunks(); ++k)
			total_samps += dedispersion_plan.get_t_processed()[i][k];
		(*output_buffer)[i] = (float **) malloc(dedispersion_plan.get_ndms()[i] * sizeof(float *));
		for (int j = 0; j < dedispersion_plan.get_ndms()[i]; ++j)
		{
			(*output_buffer)[i][j] = (float *) malloc(( total_samps ) * sizeof(float));
		}
		*output_size += ( total_samps ) * dedispersion_plan.get_ndms()[i] * sizeof(float);
	}
}

template<typename SpsParameterType>
void Sps<SpsParameterType>::allocate_memory_gpu(DedispersionPlan const &dedispersion_plan)
{
	int time_samps = dedispersion_plan.get_t_processed()[0][0] + dedispersion_plan.get_maxshift();
	_gpu_input_size = time_samps * (size_t) dedispersion_plan.get_nchans() * sizeof(unsigned short);

	printf("\n gpu inputsize: %d\n", (int)(_gpu_input_size/1024/1024));
	cudaError_t rc1 = ( cudaMalloc((void **)&_d_input, _gpu_input_size) );
	if (rc1 != cudaSuccess)
	    printf("Could not allocate memory: %d", rc1);

	if (dedispersion_plan.get_nchans() < dedispersion_plan.get_max_ndms())
	{
		_gpu_output_size = time_samps * dedispersion_plan.get_max_ndms() * sizeof(float);
	}
	else
	{
		_gpu_output_size = time_samps * dedispersion_plan.get_nchans() * sizeof(float);
	}
	printf("\n gpu outputsize: %d\n", (int)(_gpu_output_size /1024/1024));
	cudaError_t rc2 = ( cudaMalloc((void **)&_d_output, _gpu_output_size) );
	if (rc2 != cudaSuccess)
	    printf("Could not allocate memory: %d", rc2);


	( cudaMemset(_d_output, 0, _gpu_output_size) );
}

template<typename SpsParameterType>
void Sps<SpsParameterType>::operator()( unsigned device_id,
										InputData &input_data,
										DedispersionPlan &dedispersion_plan,
										UserInput &user_input,
										size_t gpu_memory,
										std::vector<float> &output_sps,
										float ****output_buffer,
										size_t *output_size)
{

	//
	long int inc = 0;
	float tstart_local = 0.0f;

	// allocate memory cpu output
	allocate_memory_cpu_output(dedispersion_plan, output_buffer, output_size);
	// allocate memory gpu
	allocate_memory_gpu(dedispersion_plan);
	//
	output_sps.resize(*output_size/sizeof(float));

	printf("\nDe-dispersing...\n");

	int t, dm_range;
	float tsamp_original = _tsamp;
	int maxshift_original = _maxshift;

	//
	float *out_tmp = NULL;
	out_tmp = (float *) malloc(( _t_processed[0][0] + _maxshift ) * _max_ndms * sizeof(float));
	memset(out_tmp, 0.0f, _t_processed[0][0] + _maxshift * _max_ndms * sizeof(float));

	// count the number of times that analysis is called - used to fill the output buffer of analysis
	int analysis_call = 0;
	for (t = 0; t < _num_tchunks; ++t)
	{

		printf("\nt_processed:\t%d, %d", _t_processed[0][t], t);

		load_data(-1, _in_bin, _d_input, &_input_buffer[(long int) ( inc * _nchans )],
							_t_processed[0][t], _maxshift, _nchans, _dmshifts);

		if (user_input.get_enable_zero_dm())
			zero_dm(_d_input, _nchans, _t_processed[0][t]+_maxshift);

		corner_turn(_d_input, _d_output, _nchans, _t_processed[0][t] + _maxshift);
		int oldBin = 1;

		for (dm_range = 0; dm_range < _range; ++dm_range)
		{
			printf("\n\n%f\t%f\t%f\t%d", _dm_low[dm_range], _dm_high[dm_range],	_dm_step[dm_range], _ndms[dm_range]), fflush(stdout);
			printf("\nAmount of telescope time processed: %f", tstart_local);

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

			_gpu_output_size = _ndms[dm_range] * ( _t_processed[dm_range][t] ) * sizeof(float);

			save_data(_d_output, out_tmp, _gpu_output_size);

			for (int k = 0; k < _ndms[dm_range]; ++k)
			{
				memcpy(&(*output_buffer)[dm_range][k][inc / _in_bin[dm_range]], &out_tmp[k * _t_processed[dm_range][t]],
							sizeof(float) * _t_processed[dm_range][t]);
			}

			if (user_input.get_output_dmt() == 1)
				write_output(dm_range, _t_processed[dm_range][t], _ndms[dm_range], gpu_memory, out_tmp, _gpu_output_size, _dm_low, _dm_high);

			if (user_input.get_enable_analysis() == 1)
			{
				analysis(dm_range, tstart_local, _t_processed[dm_range][t],
						( _t_processed[dm_range][t] + _maxshift ), _nchans,
						_maxshift, _max_ndms, _ndms, _out_bin, _sigma_cutoff,
						_d_output, _dm_low, _dm_high, _dm_step, _tsamp, output_sps, analysis_call);
				analysis_call += 1;

			}
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
