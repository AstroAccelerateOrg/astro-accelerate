#include "../Sps.h"

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
void Sps<SpsParameterType>::operator()( unsigned device_id, DataInputType const&,
                                  DataOutputType&, DedispersionPlan const &dedispersion_plan,
                                  SpsHandler, DmHandler)
{
	/*
		// Initialise the GPU.
		size_t gpu_memory=0;
		init_gpu(device_id, &gpu_memory);

		// Allocate memory on host and device.
		allocate_memory_cpu_output(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output,
	                        &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);

		// Allocate memory on host and device.
		allocate_memory_gpu(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output,
	                        &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);


		//printf("\nDe-dispersing...");
		int t, dm_range;
		double start_t, end_t;
		//start_t = omp_get_wtime();


		tsamp_original = tsamp;
		maxshift_original = maxshift;

		float *out_tmp;
		out_tmp = (float *) malloc(( t_processed[0][0] + maxshift ) * max_ndms * sizeof(float));
		memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));

		for (t = 0; t < dedispersion_plan.get_num_tchunks(); ++t)
		{
			//printf("\nt_processed:\t%d, %d", t_processed[0][t], t);
			//rfi((t_processed[0][t]+maxshift), nchans, &tmp);

			load_data(-1, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[0][t], maxshift, nchans, dmshifts);
			if (enable_zero_dm)
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

				dedisperse(dm_range, t_processed[dm_range][t], inBin, dmshifts, d_input, d_output, nchans, ( t_processed[dm_range][t] + maxshift ), maxshift, &tsamp, dm_low, dm_high, dm_step, ndms);

				gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
				//cudaDeviceSynchronize();

				save_data(d_output, out_tmp, gpu_outputsize);
				//	save_data(d_output, &output_buffer[dm_range][0][((long int)inc)/inBin[dm_range]], gpu_outputsize);

	//#pragma omp parallel for
				for (int k = 0; k < ndms[dm_range]; ++k)
				{
					memcpy(&output_buffer[dm_range][k][inc / inBin[dm_range]], &out_tmp[k * t_processed[dm_range][t]], sizeof(float) * t_processed[dm_range][t]);
				}

				if (output_dmt == 1)
					write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, out_tmp, gpu_outputsize, dm_low, dm_high);
				if (enable_analysis == 1)
					analysis(dm_range, tstart_local, t_processed[dm_range][t], ( t_processed[dm_range][t] + maxshift ), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, d_output, dm_low, dm_high, dm_step, tsamp);
				oldBin = inBin[dm_range];
			}

			memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));

			inc = inc + t_processed[0][t];
			printf("\nINC:\t%ld", inc);
			tstart_local = ( tsamp_original * inc );
			tsamp = tsamp_original;
			maxshift = maxshift_original;
		}

		cudaFree(d_input);
		cudaFree(d_output);
		free(out_tmp);
		free(input_buffer);
		free(output_buffer);*/
}

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
