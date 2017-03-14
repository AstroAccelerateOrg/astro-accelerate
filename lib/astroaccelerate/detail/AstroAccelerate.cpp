#include "../AstroAccelerate.h"

namespace ska {
namespace astroaccelerate {

template<typename AstroAccelerateParameterType>
AstroAccelerate<AstroAccelerateParameterType>::AstroAccelerate(DedispersionStrategy  &dedispersion_strategy)
						   : _num_tchunks(dedispersion_strategy.get_num_tchunks())
						   ,_range(dedispersion_strategy.get_range())
						   ,_nchans(dedispersion_strategy.get_nchans())
						   ,_maxshift(dedispersion_strategy.get_maxshift())
						   ,_tsamp(dedispersion_strategy.get_tsamp())
						   ,_max_ndms(dedispersion_strategy.get_max_ndms())
						   ,_sigma_cutoff(dedispersion_strategy.get_sigma_cutoff())
						   ,_ndms(dedispersion_strategy.get_ndms())
						   ,_dmshifts(dedispersion_strategy.get_dmshifts())
						   ,_t_processed(dedispersion_strategy.get_t_processed())
						   ,_dm_low(dedispersion_strategy.get_dm_low())
						   ,_dm_high(dedispersion_strategy.get_dm_high())
						   ,_dm_step(dedispersion_strategy.get_dm_step())
						   ,_in_bin(dedispersion_strategy.get_in_bin())
						   ,_out_bin(dedispersion_strategy.get_out_bin())
{
	_gpu_input_size = 0;
	_d_input = nullptr;
	_gpu_output_size = 0;
	_d_output = nullptr;
}

template<typename AstroAccelerateParameterType>
AstroAccelerate<AstroAccelerateParameterType>::~AstroAccelerate()
{
}

template<typename AstroAccelerateParameterType>
void AstroAccelerate<AstroAccelerateParameterType>::allocate_memory_gpu(DedispersionStrategy const &dedispersion_strategy)
{
	int time_samps = dedispersion_strategy.get_t_processed()[0][0] + dedispersion_strategy.get_maxshift();
	_gpu_input_size = time_samps * (size_t) dedispersion_strategy.get_nchans() * sizeof(unsigned short);

	//printf("\n gpu inputsize: %d\n", (int)(_gpu_input_size/1024/1024));
	cudaError_t rc1 = ( cudaMalloc((void **)&_d_input, _gpu_input_size) );
	if (rc1 != cudaSuccess)
	{
	    printf("Could not allocate gpu memory: %d", rc1);
	    exit(0);
	}

	if (dedispersion_strategy.get_nchans() < dedispersion_strategy.get_max_ndms())
	{
		_gpu_output_size = time_samps * dedispersion_strategy.get_max_ndms() * sizeof(float);
	}
	else
	{
		_gpu_output_size = time_samps * dedispersion_strategy.get_nchans() * sizeof(float);
	}
	//printf("\n gpu outputsize: %d\n", (int)(_gpu_output_size /1024/1024));
	cudaError_t rc2 = ( cudaMalloc((void **)&_d_output, _gpu_output_size) );
	if (rc2 != cudaSuccess)
	{
	    printf("Could not allocate gpu memory: %d", rc2);
	    exit(0);
	}


	( cudaMemset(_d_output, 0, _gpu_output_size) );
}

template<typename AstroAccelerateParameterType>
void AstroAccelerate<AstroAccelerateParameterType>::operator()( unsigned device_id
										,DedispersionStrategy &dedispersion_strategy
										,size_t gpu_memory
										,unsigned short *input_buffer
										,DmTime<float> &output_buffer
										,std::vector<float> &output_sps
										)
{
	//
	long int inc = 0;
	float tstart_local = 0.0f;

	// allocate memory gpu
	allocate_memory_gpu(dedispersion_strategy);
	//

	output_sps.resize(output_buffer.output_size()/sizeof(float));

	//printf("\nDe-dispersing...\n");
	GpuTimer timer;
	timer.Start();


	float tsamp_original = _tsamp;
	int maxshift_original = _maxshift;

	//float *out_tmp;
	//out_tmp = (float *) malloc(( _t_processed[0][0] + _maxshift ) * _max_ndms * sizeof(float));
	//memset(out_tmp, 0.0f, _t_processed[0][0] + _maxshift * _max_ndms * sizeof(float));

	for (int t = 0; t < _num_tchunks; ++t)
	{
		printf("\nt_processed:\t%d, %d", _t_processed[0][t], t);

		load_data(-1, _in_bin, _d_input, &input_buffer[(long int) ( inc * _nchans )], _t_processed[0][t], _maxshift, _nchans, _dmshifts);


		if (dedispersion_strategy.get_enable_zero_dm())
			zero_dm(_d_input, _nchans, _t_processed[0][t]+_maxshift);

		if (dedispersion_strategy.get_enable_zero_dm_with_outliers())
			zero_dm_outliers(_d_input, _nchans, _t_processed[0][t]+_maxshift);

		corner_turn(_d_input, _d_output, _nchans, _t_processed[0][t] + _maxshift);

		if (dedispersion_strategy.get_enable_rfi())
	 		rfi_gpu(_d_input, _nchans, _t_processed[0][t]+_maxshift);

		int oldBin = 1;
		for (int dm_range = 0; dm_range < _range; ++dm_range)
		{
			printf("\n\n%f\t%f\t%f\t%d", _dm_low[dm_range], _dm_high[dm_range], _dm_step[dm_range], _ndms[dm_range]), fflush(stdout);
			printf("\nAmount of telescope time processed: %f", tstart_local);
			_maxshift = maxshift_original / _in_bin[dm_range];

			cudaDeviceSynchronize();
			load_data(dm_range, _in_bin, _d_input, &input_buffer[(long int) ( inc * _nchans )], _t_processed[dm_range][t], _maxshift, _nchans, _dmshifts);

			if (_in_bin[dm_range] > oldBin)
			{
				bin_gpu(_d_input, _d_output, _nchans, _t_processed[dm_range - 1][t] + _maxshift * _in_bin[dm_range]);
				_tsamp = _tsamp * 2.0f;
			}

			dedisperse(dm_range, _t_processed[dm_range][t], _in_bin, _dmshifts, _d_input, _d_output, _nchans,
				( _t_processed[dm_range][t] + _maxshift ), _maxshift, &_tsamp, _dm_low, _dm_high, _dm_step, _ndms);

			if (dedispersion_strategy.get_enable_acceleration() == 1)
			{
				// gpu_outputsize = ndms[dm_range] * ( t_processed[dm_range][t] ) * sizeof(float);
				//save_data(d_output, out_tmp, gpu_outputsize);

				//#pragma omp parallel for
				for (int k = 0; k < _ndms[dm_range]; ++k)
				{
					//memcpy(&output_buffer[dm_range][k][inc / inBin[dm_range]], &out_tmp[k * t_processed[dm_range][t]], sizeof(float) * t_processed[dm_range][t]);
					save_data_offset(_d_output, k * _t_processed[dm_range][t], output_buffer[dm_range][k], inc / _in_bin[dm_range], sizeof(float) * _t_processed[dm_range][t]);
				}
				//	save_data(d_output, &output_buffer[dm_range][0][((long int)inc)/inBin[dm_range]], gpu_outputsize);
			}

//			if (output_dmt == 1)
//			{
				//for (int k = 0; k < ndms[dm_range]; k++)
				//	write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, output_buffer[dm_range][k], gpu_outputsize, dm_low, dm_high);
				//write_output(dm_range, t_processed[dm_range][t], ndms[dm_range], gpu_memory, out_tmp, gpu_outputsize, dm_low, dm_high);
//			}
			if (dedispersion_strategy.get_enable_analysis() == 1)
			{
				// TODO: put the file export back to analysis I leaving it here at the moment since for interface we need to output from the analysis.
				float *h_output_list;
				float *h_peak_list;
				size_t max_list_size, max_peak_size;
				size_t list_pos, peak_pos;
				max_list_size = (size_t) ( _ndms[dm_range]*_t_processed[dm_range][t]/2 ); // we can store 1/2 of the input plane
				max_peak_size = (size_t) ( _ndms[dm_range]*_t_processed[dm_range][t]/2 );
				h_output_list = (float*) malloc(max_list_size*4*sizeof(float)); // Allocations
				h_peak_list   = (float*) malloc(max_list_size*4*sizeof(float));

				list_pos=0;
				peak_pos=0;

				analysis_GPU(h_output_list, &list_pos, max_list_size, h_peak_list, &peak_pos, max_peak_size, dm_range, tstart_local,
							_t_processed[dm_range][t], _in_bin[dm_range], _out_bin[dm_range], &_maxshift, _max_ndms, _ndms,
							_sigma_cutoff, dedispersion_strategy.get_sigma_constant(), dedispersion_strategy.get_max_boxcar_width_in_sec()
							, _d_output, _dm_low, _dm_high, _dm_step, _tsamp);


				printf("-------> list_pos:%zu; \n", list_pos);
				//#pragma omp parallel for
				for (int count = 0; count < list_pos; count++)
				{
					h_output_list[4*count]     = h_output_list[4*count]*_dm_step[dm_range] + _dm_low[dm_range];
					h_output_list[4*count + 1] = h_output_list[4*count + 1]*_tsamp + tstart_local;
					//h_output_list[4*count + 2] = h_output_list[4*count + 2];
					//h_output_list[4*count + 3] = h_output_list[4*count + 3];
				}
				//#pragma omp parallel for
				for (int count = 0; count < peak_pos; count++)
				{
					h_peak_list[4*count]     = h_peak_list[4*count]*_dm_step[dm_range] + _dm_low[dm_range];
					h_peak_list[4*count + 1] = h_peak_list[4*count + 1]*_tsamp + tstart_local;
					//h_output_list[4*count + 2] = h_output_list[4*count + 2];
					//h_output_list[4*count + 3] = h_output_list[4*count + 3];
				}

				FILE *fp_out;
				char filename[200];

				if(list_pos>0)
				{
					sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", tstart_local, _dm_low[dm_range], _dm_high[dm_range]);
					//if ((fp_out=fopen(filename, "w")) == NULL) {
					if (( fp_out = fopen(filename, "wb") ) == nullptr)
					{
						fprintf(stderr, "Error opening output file!\n");
						exit(0);
					}
					fwrite(h_output_list, list_pos*sizeof(float), 4, fp_out);
					fclose(fp_out);
				}

				if (peak_pos > 0)
				{
					sprintf(filename, "peak_analysed-t_%.2f-dm_%.2f-%.2f.dat",tstart_local, _dm_low[dm_range],_dm_high[dm_range]);
					//if ((fp_out=fopen(filename, "w")) == NULL) {
					if ((fp_out = fopen(filename, "wb")) == nullptr)
					{
						fprintf(stderr, "Error opening output file!\n");
						exit(0);
					}
					fwrite(h_peak_list, peak_pos*sizeof(float), 4, fp_out);
					fclose(fp_out);
				}

				free(h_peak_list);
				free(h_output_list);


				// This is for testing purposes and should be removed or commented out
				//analysis_CPU(dm_range, tstart_local, t_processed[dm_range][t], (t_processed[dm_range][t]+maxshift), nchans, maxshift, max_ndms, ndms, outBin, sigma_cutoff, out_tmp,dm_low, dm_high, dm_step, tsamp);
			}
			oldBin = _in_bin[dm_range];
		}

			//memset(out_tmp, 0.0f, t_processed[0][0] + maxshift * max_ndms * sizeof(float));

		inc = inc + _t_processed[0][t];
		printf("\nINC:\t%ld", inc);
		tstart_local = ( tsamp_original * inc );
		_tsamp = tsamp_original;
		_maxshift = maxshift_original;
	}

	timer.Stop();
	float time = timer.Elapsed() / 1000;

	printf("\n\n === OVERALL DEDISPERSION THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

	printf("\n(Performed Brute-Force Dedispersion: %g (GPU estimate)",  time);
	printf("\nAmount of telescope time processed: %f", tstart_local);
	printf("\nNumber of samples processed: %ld", inc);
	printf("\nReal-time speedup factor: %lf", ( tstart_local ) / time);

	cudaFree(_d_input);
	cudaFree(_d_output);
	//free(out_tmp);
	//free(input_buffer);

	double time_processed = ( tstart_local ) / tsamp_original;
	double dm_t_processed = time_processed * dedispersion_strategy.get_total_ndms();
	double all_processed = dm_t_processed * _nchans;
	printf("\nGops based on %.2lf ops per channel per tsamp: %f", NOPS, ( ( NOPS * all_processed ) / ( time ) ) / 1000000000.0);
	int num_reg = SNUMREG;
	float num_threads = dedispersion_strategy.get_total_ndms() * ( _t_processed[0][0] ) / ( num_reg );
	float data_size_loaded = ( num_threads * _nchans * sizeof(ushort) ) / 1000000000;
	float time_in_sec = time;
	float bandwidth = data_size_loaded / time_in_sec;
	printf("\nDevice global memory bandwidth in GB/s: %f", bandwidth);
	printf("\nDevice shared memory bandwidth in GB/s: %f", bandwidth * ( num_reg ));
	float size_gb = ( _nchans * ( _t_processed[0][0] ) * sizeof(float) * 8 ) / 1000000000.0;
	printf("\nTelescope data throughput in Gb/s: %f", size_gb / time_in_sec);

	if (dedispersion_strategy.get_enable_periodicity() == 1)
	{
		//
		GpuTimer timer;
		timer.Start();
		//
		periodicity(_range, dedispersion_strategy.get_nsamp(), _max_ndms, inc, dedispersion_strategy.get_nboots(), dedispersion_strategy.get_ntrial_bins(),
					dedispersion_strategy.get_navdms(), dedispersion_strategy.get_narrow(), dedispersion_strategy.get_wide(), dedispersion_strategy.get_nsearch(),
					dedispersion_strategy.get_aggression(), dedispersion_strategy.get_sigma_cutoff(), output_buffer, _ndms, _in_bin, _dm_low, _dm_high, _dm_step, tsamp_original);
		//
		timer.Stop();
		float time = timer.Elapsed()/1000;
		printf("\n\n === OVERALL PERIODICITY THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Peroidicity Location: %f (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", inc);
		printf("\nReal-time speedup factor: %f", ( tstart_local ) / ( time ));
	}

	if (dedispersion_strategy.get_enable_acceleration() == 1)
	{
		// Input needed for fdas is output_buffer which is DDPlan
		// Assumption: gpu memory is free and available
		//
		GpuTimer timer;
		timer.Start();
		// acceleration(range, nsamp, max_ndms, inc, nboots, ntrial_bins, navdms, narrow, wide, nsearch, aggression, sigma_cutoff, output_buffer, ndms, inBin, dm_low, dm_high, dm_step, tsamp_original);
		acceleration_fdas(_range
						  ,dedispersion_strategy.get_nsamp()
						  ,_max_ndms
						  ,inc
						  ,dedispersion_strategy.get_nboots()
						  ,dedispersion_strategy.get_ntrial_bins()
						  ,dedispersion_strategy.get_navdms()
						  ,dedispersion_strategy.get_narrow()
						  ,dedispersion_strategy.get_wide()
						  ,dedispersion_strategy.get_nsearch()
						  ,dedispersion_strategy.get_aggression()
						  ,dedispersion_strategy.get_sigma_cutoff()
						  ,output_buffer
						  ,_ndms
						  ,_in_bin
						  ,_dm_low
						  ,_dm_high
						  ,_dm_step
						  ,tsamp_original
						  ,dedispersion_strategy.get_enable_fdas_custom_fft()
						  ,dedispersion_strategy.get_enable_fdas_inbin()
						  ,dedispersion_strategy.get_enable_fdas_norm());
		//
		timer.Stop();
		float time = timer.Elapsed()/1000;
		printf("\n\n === OVERALL TDAS THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

		printf("\nPerformed Acceleration Location: %lf (GPU estimate)", time);
		printf("\nAmount of telescope time processed: %f", tstart_local);
		printf("\nNumber of samples processed: %ld", inc);
		printf("\nReal-time speedup factor: %lf", ( tstart_local ) / ( time ));
	}
	//free(out_tmp);
	cudaFree(_d_input);
	cudaFree(_d_output);
}

} // namespace astroaccelerate
}
