#include "../Sps.h"
#include <algorithm>
#include <vector>

namespace ska {
namespace astroaccelerate {
namespace sps {

template<typename SpsParameterType>
Sps<SpsParameterType>::Sps()
{
}

template<typename SpsParameterType>
Sps<SpsParameterType>::~Sps()
{
}

template<typename SpsParameterType>
template<typename SpsHandler, typename DmHandler>
void Sps<SpsParameterType>::operator()( unsigned device_id, IOData &io_data, DedispersionPlan &dedispersion_plan,
                                        UserInput const &user_input, SpsHandler, DmHandler)
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

		// Allocate memory on host and device.
		io_data.allocate_memory_cpu_output(dedispersion_plan);
		io_data.allocate_memory_gpu(dedispersion_plan);

		//printf("\nDe-dispersing...");
		int t, dm_range;
		//double start_t, end_t;
		//start_t = omp_get_wtime();

		//
		int tsamp_original = dedispersion_plan.get_tsamp();
		int maxshift = dedispersion_plan.get_maxshift();
		int max_ndms = dedispersion_plan.get_max_ndms();
		int maxshift_original = maxshift;
		auto t_processed = dedispersion_plan.get_t_processed();
		int num_tchunks = dedispersion_plan.get_num_tchunks();
		int nchans = dedispersion_plan.get_nchans();
		float* dmshifts = dedispersion_plan.get_dmshifts();
		unsigned short* input_buffer = io_data.get_input_buffer();
		auto output_buffer = io_data.get_output_buffer();
		int* inBin = dedispersion_plan.get_in_bin();
		int* outBin = dedispersion_plan.get_out_bin();
		unsigned short* d_input = io_data.get_d_input();
		float* d_output = io_data.get_d_output();
		int range = dedispersion_plan.get_range();
		float tsamp = dedispersion_plan.get_tsamp();
		float* dm_low  = dedispersion_plan.get_dm_low();
		float* dm_high = dedispersion_plan.get_dm_high();
		float* dm_step = dedispersion_plan.get_dm_step();
		int* ndms = dedispersion_plan.get_ndms();
		size_t gpu_outputsize = io_data.get_gpu_output_size();
		float sigma_cutoff = user_input.get_sigma_cutoff();

		//
		std::vector<float> out_tmp(t_processed[0][0] + maxshift * max_ndms, 0.0f);

		for (t = 0; t < num_tchunks; ++t)
		{
			//printf("\nt_processed:\t%d, %d", t_processed[0][t], t);
			//rfi((t_processed[0][t]+maxshift), nchans, &tmp);


			load_data(-1, inBin, d_input, &input_buffer[(long int) ( inc * nchans )],
								t_processed[0][t], maxshift, nchans, dmshifts);

			if (user_input.get_enable_zero_dm())
				zero_dm(d_input, nchans, t_processed[0][t]+maxshift);

			corner_turn(d_input, d_output, nchans, t_processed[0][t] + maxshift);
			int oldBin = 1;

			for (dm_range = 0; dm_range < range; ++dm_range)
			{

				//printf("\n\n%f\t%f\t%f\t%d", dm_low[dm_range], dm_high[dm_range], dm_step[dm_range], ndms[dm_range]), fflush(stdout);
				//printf("\nAmount of telescope time processed: %f", tstart_local);
				maxshift = maxshift_original / inBin[dm_range];

				cudaDeviceSynchronize();
				load_data(dm_range, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[dm_range][t], maxshift, nchans, dmshifts);

				if (inBin[dm_range] > oldBin)
				{
					bin_gpu(d_input, d_output, nchans, t_processed[dm_range - 1][t] + maxshift * inBin[dm_range]);
					( tsamp ) = ( tsamp ) * 2.0f;
				}

				dedisperse(dm_range, t_processed[dm_range][t], inBin, dmshifts, d_input, d_output, nchans,
				           ( t_processed[dm_range][t] + maxshift ), maxshift, &tsamp, dm_low, dm_high, dm_step, ndms);

				gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
				//cudaDeviceSynchronize();

				save_data(d_output, out_tmp.data(), gpu_outputsize);
				//	save_data(d_output, &output_buffer[dm_range][0][((long int)inc)/inBin[dm_range]], gpu_outputsize);


				//#pragma omp parallel for
				for (int k = 0; k < ndms[dm_range]; ++k)
				{
					memcpy(&output_buffer[dm_range][k][inc / inBin[dm_range]], &out_tmp[k * t_processed[dm_range][t]],
								sizeof(float) * t_processed[dm_range][t]);
				}

				if (user_input.get_output_dmt() == 1)
					write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, out_tmp.data(), gpu_outputsize, dm_low, dm_high);
				if (user_input.get_enable_analysis() == 1)
					analysis(dm_range, tstart_local, t_processed[dm_range][t], ( t_processed[dm_range][t] + maxshift ), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, d_output, dm_low, dm_high, dm_step, tsamp);
				oldBin = inBin[dm_range];
			}

			std::fill(out_tmp.begin(), out_tmp.end(), 0.0f);

			inc = inc + t_processed[0][t];
			printf("\nINC:\t%ld", inc);
			tstart_local = ( tsamp_original * inc );
			tsamp = tsamp_original;
			maxshift = maxshift_original;
		}

		cudaFree(d_input);
		cudaFree(d_output);
}

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
