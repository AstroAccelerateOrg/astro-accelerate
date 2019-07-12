#ifndef ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_GENERIC_HPP
#define ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_GENERIC_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>
#include "aa_pipeline.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_ddtr_plan.hpp"
#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"
#include "aa_periodicity_strategy.hpp"
#include "aa_fdas_strategy.hpp"

#include "aa_filterbank_metadata.hpp"
#include "aa_device_load_data.hpp"
#include "aa_bin_gpu.hpp"
#include "aa_zero_dm.hpp"
#include "aa_zero_dm_outliers.hpp"
#include "aa_corner_turn.hpp"
#include "aa_device_rfi.hpp"
#include "aa_dedisperse.hpp"

#include "aa_gpu_timer.hpp"

#include "aa_device_analysis.hpp"
#include "aa_device_periods.hpp"
#include "aa_device_acceleration_fdas.hpp"
#include "aa_pipeline_runner.hpp"

#include "aa_gpu_timer.hpp"

namespace astroaccelerate {

	/**
	* \class aa_permitted_pipelines_5 aa_permitted_pipelines_5.hpp "include/aa_permitted_pipelines_5.hpp"
	* \brief Templated class to run dedispersion and analysis and periodicity and acceleration.
	* \details The class is templated over the zero_dm_type (aa_pipeline::component_option::zero_dm or aa_pipeline::component_option::zero_dm_with_outliers).
	* \author AstroAccelerate team.
	*/

	//template<aa_pipeline::component_option zero_dm_type, bool enable_old_rfi>
	class aa_permitted_pipelines_generic : public aa_pipeline_runner {
	private:
		float              ***m_output_buffer;
		int                **t_processed; //This should be in ddtr_strategy
		const aa_pipeline::pipeline          m_pipeline_components; /** The user requested pipeline that was bound to the aa_pipeline_api instance on construction. */
		const aa_pipeline::pipeline_option   m_pipeline_options; /** The user requested pipeline details containing component options for the aa_pipeline_api instance. */
		aa_ddtr_strategy                     m_ddtr_strategy;
		aa_analysis_strategy                 m_analysis_strategy;
		aa_periodicity_strategy              m_periodicity_strategy;
		aa_fdas_strategy                     m_fdas_strategy;
		unsigned short     const*const m_input_buffer;
		int                num_tchunks; 
		std::vector<float> dm_shifts; 
		float              *dmshifts; 
		int                maxshift; 
		int                max_ndms; 
		int                nchans; 
		int                nbits;
		int                failsafe;
		long int           inc;
		float              tsamp;
		float              tsamp_original;
		int                maxshift_original;
		int                nRanges;
		float              tstart_local;
		int                oldBin;

		unsigned short     *d_DDTR_input;
		float              *d_DDTR_output;

		std::vector<float> dm_low;
		std::vector<float> dm_high;
		std::vector<float> dm_step;
		std::vector<int>   inBin;
		
		// single pulse detection output candidates
		unsigned int *h_SPD_candidate_list_DM;
		unsigned int *h_SPD_candidate_list_TS;
		float        *h_SPD_candidate_list_SNR;
		unsigned int *h_SPD_candidate_list_BW;
		size_t        SPD_max_peak_size;
		size_t        SPD_nCandidates;
		

		// fdas acceleration search settings
		bool m_fdas_enable_custom_fft;
		bool m_fdas_enable_inbin;
		bool m_fdas_enable_norm;
		bool m_fdas_enable_output_ffdot_plan;
		bool m_fdas_enable_output_list;

		// Pipeline flags
		bool do_dedispersion;
		bool do_single_pulse_detection;
		bool do_periodicity_search;
		bool do_fdas;
		
		bool do_copy_DDTR_data_to_host;
		
		bool memory_allocated;
		bool memory_cleanup;
		bool periodicity_did_run;
		bool acceleration_did_run;
		bool did_notify_of_finishing_component;

		//Loop counter variables
		int current_time_chunk;
		int current_range;
		aa_gpu_timer       m_timer;

		float  *m_d_MSD_workarea = NULL;
		float  *m_d_MSD_interpolated = NULL;
		ushort *m_d_MSD_output_taps = NULL;

		bool cleanup_DDTR() {
			if (memory_allocated && !memory_cleanup) {
				LOG(log_level::debug, "DDTR -> Memory cleanup after de-dispersion");
				cudaFree(d_DDTR_input);
				cudaFree(d_DDTR_output);
				
				if(do_single_pulse_detection){
					cudaFree(m_d_MSD_workarea);
					cudaFree(m_d_MSD_output_taps);
					cudaFree(m_d_MSD_interpolated);
					
					free(h_SPD_candidate_list_DM);
					free(h_SPD_candidate_list_TS);
					free(h_SPD_candidate_list_SNR);
					free(h_SPD_candidate_list_BW);
				}
				
				// Why this is not in the ddtr_strategy?
				size_t t_processed_size = m_ddtr_strategy.t_processed().size();
				for (size_t i = 0; i < t_processed_size; i++) {
					free(t_processed[i]);
				}
				free(t_processed);

				memory_cleanup = false;
			}
			return true;
		}
		
		/** \brief Allocate the GPU memory needed for dedispersion. */
		void allocate_gpu_memory_DDTR(const int &maxshift, const int &max_ndms, const int &nchans, int **const t_processed, unsigned short **const d_DDTR_input, float **const d_DDTR_output) {
			int time_samps = t_processed[0][0] + maxshift;
			printf("\n\n\n%d\n\n\n", time_samps);
			size_t gpu_inputsize = (size_t)time_samps * (size_t)nchans * sizeof(unsigned short);
			cudaError_t e;

			e = cudaMalloc((void **)d_DDTR_input, gpu_inputsize);
			if (e != cudaSuccess) {
				LOG(log_level::error, "Could not allocate memory for d_DDTR_input using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}

			size_t gpu_outputsize = 0;
			if (nchans < max_ndms) {
				gpu_outputsize = (size_t)time_samps * (size_t)max_ndms * sizeof(float);
			}
			else {
				gpu_outputsize = (size_t)time_samps * (size_t)nchans * sizeof(float);
			}

			e = cudaMalloc((void **)d_DDTR_output, gpu_outputsize);
			if (e != cudaSuccess) {
				LOG(log_level::error, "Could not allocate memory for d_DDTR_output using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}

			cudaMemset(*d_DDTR_output, 0, gpu_outputsize);
		}


		/** \brief Allocate memory for single pulse detection (SPD) on the host. */
		void allocate_cpu_memory_SPD(size_t max_nDMs, size_t timesamples_per_chunk){
			SPD_max_peak_size        = (size_t)(max_nDMs*timesamples_per_chunk/2);
			
			h_SPD_candidate_list_DM  = (unsigned int*)malloc(SPD_max_peak_size*sizeof(unsigned int));
			h_SPD_candidate_list_TS  = (unsigned int*)malloc(SPD_max_peak_size*sizeof(unsigned int));
			h_SPD_candidate_list_SNR = (float*)malloc(SPD_max_peak_size*sizeof(float));
			h_SPD_candidate_list_BW  = (unsigned int*)malloc(SPD_max_peak_size*sizeof(unsigned int));
			if(h_SPD_candidate_list_DM==NULL || h_SPD_candidate_list_TS==NULL || h_SPD_candidate_list_SNR==NULL || h_SPD_candidate_list_BW==NULL) {
				LOG(log_level::error, "Could not allocate memory on the host for single pulse detection candidates");
			}
		}
		
		/** \brief Allocate memory for single pulse detection (SPD) on the device.	*/
		void allocate_gpu_memory_SPD(float **const d_MSD_workarea, unsigned short **const d_MSD_output_taps, float **const d_MSD_interpolated, const unsigned long int &MSD_maxtimesamples, const size_t &MSD_profile_size) {
			cudaError_t e;
			
			e = cudaMalloc((void **)d_MSD_workarea, MSD_maxtimesamples*5.5*sizeof(float));
			if (e != cudaSuccess) {
				LOG(log_level::error, "Could not allocate memory for d_MSD_workarea using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}

			e = cudaMalloc((void **) &(*d_MSD_output_taps), sizeof(ushort)*2*MSD_maxtimesamples);
			if (e != cudaSuccess) {
				LOG(log_level::error, "Could not allocate memory for d_MSD_output_taps using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}

			e = cudaMalloc((void **)d_MSD_interpolated, sizeof(float)*MSD_profile_size);
			if (e != cudaSuccess) {
				LOG(log_level::error, "Could not allocate memory for d_MSD_interpolated cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}

		}

		/**
		* \brief Allocate a 3D array that is an output buffer that stores dedispersed array data.
		* \details This array is used by periodicity.
		*/
		void allocate_memory_cpu_output() {
			size_t estimate_outputbuffer_size = 0;
			size_t outputsize = 0;
			const size_t nRanges = m_ddtr_strategy.get_nRanges();
			const int *ndms = m_ddtr_strategy.ndms_data();

			for (size_t i = 0; i < nRanges; i++) {
				for (int j = 0; j < m_ddtr_strategy.num_tchunks(); j++) {
					estimate_outputbuffer_size += (size_t)(t_processed[i][j]*sizeof(float)*ndms[i]);
				}
			}

			//if(do_copy_DDTR_data_to_host){
				outputsize = 0;
				m_output_buffer = (float ***)malloc(nRanges * sizeof(float **));
				for (size_t i = 0; i < nRanges; i++) {
					int total_samps = 0;
					for (int k = 0; k < num_tchunks; k++) {
						total_samps += t_processed[i][k];
					}
					m_output_buffer[i] = (float **)malloc(ndms[i] * sizeof(float *));
					for (int j = 0; j < ndms[i]; j++) {
						m_output_buffer[i][j] = (float *)malloc((total_samps) * sizeof(float));
					}
					outputsize += (total_samps)* ndms[i] * sizeof(float);
				}
			//}
		}

		/** \brief Method that allocates all memory for this pipeline. */
		bool set_data() {
			num_tchunks = m_ddtr_strategy.num_tchunks();
			size_t t_processed_size = m_ddtr_strategy.t_processed().size();

			t_processed = new int*[t_processed_size];
			for (size_t i = 0; i < t_processed_size; i++) {
				t_processed[i] = new int[m_ddtr_strategy.t_processed().at(i).size()];
			}

			for (size_t i = 0; i < t_processed_size; i++) {
				for (size_t j = 0; j < m_ddtr_strategy.t_processed().at(i).size(); j++) {
					t_processed[i][j] = m_ddtr_strategy.t_processed().at(i).at(j);
				}
			}

			dm_shifts = m_ddtr_strategy.dmshifts();
			dmshifts = dm_shifts.data();
			maxshift = m_ddtr_strategy.maxshift();
			max_ndms = m_ddtr_strategy.max_ndms();
			nchans = m_ddtr_strategy.metadata().nchans();
			nbits = m_ddtr_strategy.metadata().nbits();
			failsafe = 0;
			inc = 0;
			tsamp_original = m_ddtr_strategy.metadata().tsamp();
			maxshift_original = maxshift;
			nRanges = m_ddtr_strategy.get_nRanges();
			tstart_local = 0.0;

			//Data for dedispersion
			d_DDTR_input = NULL;
			d_DDTR_output = NULL;

			//Allocate GPU memory for dedispersion
			allocate_gpu_memory_DDTR(maxshift, max_ndms, nchans, t_processed, &d_DDTR_input, &d_DDTR_output);
			
			//Allocate GPU memory for SPD (i.e. analysis)
			if(do_single_pulse_detection) {
				allocate_cpu_memory_SPD(max_ndms, t_processed[0][0]);
				allocate_gpu_memory_SPD(&m_d_MSD_workarea, &m_d_MSD_output_taps, &m_d_MSD_interpolated, m_analysis_strategy.MSD_data_info(), m_analysis_strategy.MSD_profile_size_in_bytes());
			}
			//Allocate memory for CPU output for periodicity
			allocate_memory_cpu_output();

			//Put the dm low, high, step struct contents into separate arrays again.
			//This is needed so that the kernel wrapper functions don't need to be modified.
			dm_low.resize(m_ddtr_strategy.get_nRanges());
			dm_high.resize(m_ddtr_strategy.get_nRanges());
			dm_step.resize(m_ddtr_strategy.get_nRanges());
			inBin.resize(m_ddtr_strategy.get_nRanges());
			for (size_t i = 0; i < m_ddtr_strategy.get_nRanges(); i++) {
				dm_low[i] = m_ddtr_strategy.dm(i).low;
				dm_high[i] = m_ddtr_strategy.dm(i).high;
				dm_step[i] = m_ddtr_strategy.dm(i).step;
				inBin[i] = m_ddtr_strategy.dm(i).inBin;
			}
			memory_allocated = true;
			return true;
		}

		/** \brief Method that pipeline flags according to requested pipeline components. */
		void set_pipeline_flags(){
			//----> Components
			//const aa_pipeline::component cmp_dedispersion = aa_pipeline::component::dedispersion;
			const aa_pipeline::component cmp_analysis     = aa_pipeline::component::analysis;
			const aa_pipeline::component cmp_periodicity  = aa_pipeline::component::periodicity;
			const aa_pipeline::component cmp_fdas         = aa_pipeline::component::fdas;
			
			//----> Component options
			const aa_pipeline::component_option opt_copy_ddtr_data_to_host = aa_pipeline::component_option::copy_ddtr_data_to_host;
			
			do_dedispersion = true; // default component
			if(m_pipeline_components.find(cmp_analysis) != m_pipeline_components.end()) do_single_pulse_detection = true; 
			else do_single_pulse_detection = false;
			if(m_pipeline_components.find(cmp_periodicity) != m_pipeline_components.end()) do_periodicity_search = true;
			else do_periodicity_search = false;
			if(m_pipeline_components.find(cmp_fdas) != m_pipeline_components.end()) do_fdas = true;
			else do_fdas = false;
			
			do_copy_DDTR_data_to_host = false;
			if(m_pipeline_options.find(opt_copy_ddtr_data_to_host) != m_pipeline_options.end()) do_copy_DDTR_data_to_host = true;
			if(do_periodicity_search) do_copy_DDTR_data_to_host = true;
			if(do_fdas) do_copy_DDTR_data_to_host = true;
		}

		/** \brief Transfer data from the device to the host. */
		inline void save_data_offset(float *device_pointer, int device_offset, float *host_pointer, int host_offset, size_t size) {
			cudaMemcpy(host_pointer + host_offset, device_pointer + device_offset, size, cudaMemcpyDeviceToHost);
		}

		/** \brief Transfer data from the device to the host. */
		inline void save_data(float *const device_pointer, float *const host_pointer, const size_t &size) {
			cudaMemcpy(host_pointer, device_pointer, size, cudaMemcpyDeviceToHost);
		}

		/**
		* \brief Run the pipeline by processing the next time chunk of data.
		* \details Process any flags for dumping output or providing it back to the user.
		* \details aa_permitted_pipelines_5 does not return intermediate data.
		* \details If this is required, then optional parameters can be passed as arguments which would contain periodicity and acceleration output.
		* \details Since it is the pipelines perogative to decide when periodicity and acceleration are run, the user is responsible for checking when these objects were actually modified by the pipeline.
		* \details The user is expected to keep calling run_pipeline until it returns false.
		* \details Returning true indicates there is more to process.
		* \details Returning false indicates the pipeline is finished running.
		* \returns A boolean to indicate whether further time chunks are available to process (true) or not (false).
		*/
		bool run_pipeline(aa_pipeline_runner::status &status_code) {
			const aa_pipeline::component_option opt_zero_dm                = aa_pipeline::component_option::zero_dm;
			const aa_pipeline::component_option opt_zero_dm_with_outliers  = aa_pipeline::component_option::zero_dm_with_outliers;
			const aa_pipeline::component_option opt_old_rfi                = aa_pipeline::component_option::old_rfi;

			printf("NOTICE: Pipeline start/resume run_pipeline_generic.\n");
			if (current_time_chunk >= num_tchunks) {
				
				//------------------> End of DDTR + SPD
				if (!did_notify_of_finishing_component) {
					m_timer.Stop();
					float time = m_timer.Elapsed() / 1000;
					printf("\n\n === OVERALL DEDISPERSION THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");
					printf("\n(Performed Brute-Force Dedispersion: %g (GPU estimate)", time);
					printf("\nAmount of telescope time processed: %f", tstart_local);
					printf("\nNumber of samples processed: %ld", inc);
					printf("\nReal-time speedup factor: %lf\n", (tstart_local) / time);
					
					cleanup_DDTR();
					status_code = aa_pipeline_runner::status::finished_component;
					did_notify_of_finishing_component = true;
					return true;
				}
				//--------------------------------------------------------------------------------<

				
				//------------------> Periodicity
				if (!periodicity_did_run && do_periodicity_search) {
					bool periodicity_return_value = periodicity();
					status_code = aa_pipeline_runner::status::finished_component;
					return periodicity_return_value;
				}
				//--------------------------------------------------------------------------------<

				
				//------------------> Acceleration FDAS
				if (!acceleration_did_run && do_fdas) {
					bool acceleration_return_value = acceleration();
					status_code = aa_pipeline_runner::status::finished;
					return acceleration_return_value;
				}
				//--------------------------------------------------------------------------------<

				return false; // In this case, there are no more chunks to process, and periodicity and acceleration both ran.
			}
			else if (current_time_chunk == 0) {
				m_timer.Start();
			}
			
			//----------------------------------------------------------------------->
			//------------------> Dedispersion
			const int *ndms = m_ddtr_strategy.ndms_data();
			
			if(current_range==0) {
				printf("\nNOTICE: t_processed:\t%d, %d", t_processed[0][current_time_chunk], current_time_chunk);

				//checkCudaErrors(cudaGetLastError());
				load_chunk_data(d_DDTR_input, &m_input_buffer[(long int)(inc * nchans)], t_processed[0][current_time_chunk], maxshift_original, nchans, dmshifts);
				//checkCudaErrors(cudaGetLastError());
				
				//---> Zero DM
				if (m_pipeline_options.find(opt_zero_dm) != m_pipeline_options.end()) {
					printf("\nPerforming zero DM...");
					zero_dm(d_DDTR_input, nchans, t_processed[0][current_time_chunk] + maxshift_original, nbits);
				}
				else if (m_pipeline_options.find(opt_zero_dm_with_outliers) != m_pipeline_options.end()) {
					printf("\nPerforming zero dM with outliers...");
					zero_dm_outliers(d_DDTR_input, nchans, t_processed[0][current_time_chunk] + maxshift_original);
				}
				//-------------------<

	
				//checkCudaErrors(cudaGetLastError());

				corner_turn(d_DDTR_input, d_DDTR_output, nchans, t_processed[0][current_time_chunk] + maxshift_original);

				//checkCudaErrors(cudaGetLastError());

				//---> old RFI
				if (m_pipeline_options.find(opt_old_rfi) != m_pipeline_options.end()) {
					printf("\nPerforming old GPU rfi...");
					rfi_gpu(d_DDTR_input, nchans, t_processed[0][current_time_chunk] + maxshift_original);
				}
				//-------------------<

				//checkCudaErrors(cudaGetLastError());
				oldBin = 1;
			}

			
			//for (size_t dm_range = 0; dm_range < range; dm_range++) {
			if(current_time_chunk<num_tchunks){
				int dm_range = current_range;
				float tsamp = tsamp_original*((float) inBin[dm_range]);
				printf("\n\nNOTICE: %f\t%f\t%f\t%d\n", m_ddtr_strategy.dm(dm_range).low, m_ddtr_strategy.dm(dm_range).high, m_ddtr_strategy.dm(dm_range).step, m_ddtr_strategy.ndms(dm_range));
				printf("\nAmount of telescope time processed: %f\n", tstart_local);

				maxshift = maxshift_original / inBin[dm_range];

				cudaDeviceSynchronize();
				//checkCudaErrors(cudaGetLastError());

				set_dedispersion_constants(t_processed[dm_range][current_time_chunk], maxshift);

				//checkCudaErrors(cudaGetLastError());

				
				if (inBin[dm_range] > oldBin) {
					bin_gpu(d_DDTR_input, d_DDTR_output, nchans, t_processed[dm_range - 1][current_time_chunk] + maxshift * inBin[dm_range]);
				}

				//checkCudaErrors(cudaGetLastError());
				
				dedisperse(dm_range, t_processed[dm_range][current_time_chunk], inBin.data(), dmshifts, d_DDTR_input, d_DDTR_output, nchans, &tsamp, dm_low.data(), dm_step.data(), ndms, nbits, failsafe);

				//checkCudaErrors(cudaGetLastError());
				
				//-----------> Copy data to the host
				if(do_copy_DDTR_data_to_host){
					for (int k = 0; k < ndms[dm_range]; k++) {
						save_data_offset(d_DDTR_output, k * t_processed[dm_range][current_time_chunk], m_output_buffer[dm_range][k], inc / inBin[dm_range], sizeof(float) * t_processed[dm_range][current_time_chunk]);
					}
				}
				//--------------------------------------------<
				

				//------------------> Single pulse detection
				if (do_single_pulse_detection) {
					const bool dump_to_disk = true;
					const bool dump_to_user = false;
					analysis_output output;
					SPD_nCandidates = 0;
					analysis_GPU(
						h_SPD_candidate_list_DM,
						h_SPD_candidate_list_TS,
						h_SPD_candidate_list_SNR,
						h_SPD_candidate_list_BW,
						&SPD_nCandidates,
						SPD_max_peak_size,
						dm_range,
						tstart_local,
						t_processed[dm_range][current_time_chunk],
						inBin[dm_range],
						&maxshift,
						max_ndms,
						ndms,
						m_analysis_strategy.sigma_cutoff(),
						m_analysis_strategy.sigma_constant(),
						m_analysis_strategy.max_boxcar_width_in_sec(),
						d_DDTR_output,
						dm_low.data(),
						dm_high.data(),
						dm_step.data(),
						tsamp,
						m_analysis_strategy.candidate_algorithm(),
						m_d_MSD_workarea,
						m_d_MSD_output_taps,
						m_d_MSD_interpolated,
						m_analysis_strategy.MSD_data_info(),
						m_analysis_strategy.enable_msd_baseline_noise(),
						dump_to_disk,
						dump_to_user,
						output);
				}
				//--------------------------------------------------------------------------------<
				oldBin = inBin[dm_range];
			}

			printf("NOTICE: Pipeline ended run_pipeline_5 over chunk %d / %d and range %d / %d.\n", current_time_chunk, num_tchunks, current_range, nRanges);
			
			++current_range;
			if(current_range>=nRanges) {
				inc = inc + t_processed[0][current_time_chunk];
				tstart_local = (tsamp_original * inc);
				printf("\nNOTICE: INC:\t%ld\n", inc);
				++current_time_chunk;
				current_range = 0;
			}
			
			status_code = aa_pipeline_runner::status::has_more;
			return true;
		}

		bool periodicity() {
			if (periodicity_did_run) return false;
			aa_gpu_timer timer;
			timer.Start();
			const int *ndms = m_ddtr_strategy.ndms_data();
			GPU_periodicity(m_ddtr_strategy.get_nRanges(),
				m_ddtr_strategy.metadata().nsamples(),
				max_ndms,
				inc,
				m_periodicity_strategy.sigma_cutoff(),
				m_output_buffer,
				ndms,
				inBin.data(),
				dm_low.data(),
				dm_high.data(),
				dm_step.data(),
				tsamp_original,
				m_periodicity_strategy.nHarmonics(),
				m_periodicity_strategy.candidate_algorithm(),
				m_periodicity_strategy.enable_msd_baseline_noise(),
				m_periodicity_strategy.sigma_constant());

			timer.Stop();
			float time = timer.Elapsed()/1000;
			printf("\n\n === OVERALL PERIODICITY THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

			printf("\nPerformed Periodicity Location: %f (GPU estimate)", time);
			printf("\nAmount of telescope time processed: %f", tstart_local);
			printf("\nNumber of samples processed: %ld", inc);
			printf("\nReal-time speedup factor: %f\n", (tstart_local) / (time));
			periodicity_did_run = true;
			return true;
		}

		bool acceleration() {
			// Input needed for fdas is output_buffer.
			// Assumption: GPU memory is free and available.
			if (acceleration_did_run) return false;
			aa_gpu_timer timer;
			timer.Start();
			const int *ndms = m_ddtr_strategy.ndms_data();

			acceleration_fdas(m_ddtr_strategy.get_nRanges(),
				m_ddtr_strategy.metadata().nsamples(),
				max_ndms,
				inc,
				m_fdas_strategy.sigma_cutoff(),
				m_output_buffer,
				ndms,
				inBin.data(),
				dm_low.data(),
				dm_high.data(),
				dm_step.data(),
				tsamp_original,
				m_fdas_enable_custom_fft,
				m_fdas_enable_inbin,
				m_fdas_enable_norm,
				m_fdas_strategy.sigma_constant(),
				m_fdas_enable_output_ffdot_plan,
				m_fdas_enable_output_list);

			timer.Stop();
			float time = timer.Elapsed()/1000;
			printf("\n\n === OVERALL TDAS THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

			printf("\nPerformed Acceleration Location: %lf (GPU estimate)", time);
			printf("\nAmount of telescope time processed: %f", tstart_local);
			printf("\nNumber of samples processed: %ld", inc);
			printf("\nReal-time speedup factor: %lf\n", (tstart_local) / (time));
			acceleration_did_run = true;

			/*
			for (size_t i = 0; i < range; i++) {
				for (int j = 0; j < ndms[i]; j++) {
					free(m_output_buffer[i][j]);
				}
				free(m_output_buffer[i]);
			}
			free(m_output_buffer);
			*/

			return true;
		}

	public:
		aa_permitted_pipelines_generic(
		const aa_pipeline::pipeline &pipeline_components,
		const aa_pipeline::pipeline_option &pipeline_options,
		const aa_ddtr_strategy &ddtr_strategy,
		const aa_analysis_strategy &analysis_strategy,
		const aa_periodicity_strategy &periodicity_strategy,
		const aa_fdas_strategy &fdas_strategy,
		const bool &fdas_enable_custom_fft,
		const bool &fdas_enable_inbin,
		const bool &fdas_enable_norm,
		const bool &fdas_enable_output_ffdot_plan,
		const bool &fdas_enable_output_list,
		unsigned short const*const input_buffer)
		:
		m_pipeline_components(pipeline_components),
		m_pipeline_options(pipeline_options),
		m_ddtr_strategy(ddtr_strategy),
		m_analysis_strategy(analysis_strategy),
		m_periodicity_strategy(periodicity_strategy),
		m_fdas_strategy(fdas_strategy),
		m_input_buffer(input_buffer),
		m_fdas_enable_custom_fft(fdas_enable_custom_fft),
		m_fdas_enable_inbin(fdas_enable_inbin),
		m_fdas_enable_norm(fdas_enable_norm),
		m_fdas_enable_output_ffdot_plan(fdas_enable_output_ffdot_plan),
		m_fdas_enable_output_list(fdas_enable_output_list),
		memory_allocated(false),
		memory_cleanup(false),
		periodicity_did_run(false),
		acceleration_did_run(false),
		did_notify_of_finishing_component(false),
		current_time_chunk(0),
		current_range(0),
		m_d_MSD_workarea(NULL),
		m_d_MSD_interpolated(NULL),
		m_d_MSD_output_taps(NULL) {
	}

		~aa_permitted_pipelines_generic() {
			//Only call cleanup if memory had been allocated during setup,
			//and if the memory was not already cleaned up usingthe cleanup method.
			if (memory_allocated && !memory_cleanup) {
				cleanup();
			}
		}

		aa_permitted_pipelines_generic(const aa_permitted_pipelines_generic &) = delete;
		
		void printinfo(){
			printf("----------------------------------------------------------\n");
			printf("Pipeline components:\n");
			for (auto const i : m_pipeline_components) {
				//LOG(log_level::notice, component_name(i));
				printf("%s\n", component_name(i).c_str());
			}
			printf("\n");
			printf("Component options:\n");
			for (auto const i : m_pipeline_options) {
				//LOG(log_level::notice, component_option_description(i));
				printf("%s\n", component_option_description(i).c_str());
			}
			printf("----------------------------------------------------------\n");
		}
		
		/** \brief Method to setup and allocate memory for the pipeline containers. */
		bool setup() override {
			if (!memory_allocated) {
				set_pipeline_flags();
				return set_data();
			}

			if (memory_allocated) {
				return true;
			}

			return false;
		}

		/** \brief Override base class next() method to process next time chunk. */
		bool next() override {
			if (memory_allocated) {
				aa_pipeline_runner::status tmp;
				return run_pipeline(tmp);
			}
			return false;
		}

		/** \brief Override base class next() method to process next time chunk. Also provides a status code. */
		bool next(aa_pipeline_runner::status &status_code) override {
			if (memory_allocated) {
				return run_pipeline(status_code);
			}
			return false;
		}


		/**
		 * \brief Return the pointer to the complete dedispersed output data.
		 * \details The array data is only useful once the pipeline has finished running.
		 * \details Users should finish running the pipeline so that all dedispersion output is available.
		 * \details Alternatively, the user may access the data one time chunk at a time, but the next time chunk will not have been computed yet.
		 * \details The structure of the ddtr output buffer data is indexed by: time_chunk index, dm_range index, dm.
		 */
		float*** output_buffer() {
			if(memory_allocated) {
				return m_output_buffer;
			}
			return NULL;
		}

		float* h_SPD_snr(){
			if(memory_allocated && !memory_cleanup){
				return h_SPD_candidate_list_SNR;
			}
			return NULL;
		}

		unsigned int* h_SPD_ts(){
			if(memory_allocated && !memory_cleanup){
					return h_SPD_candidate_list_TS;
			}
			return NULL;
		}

		unsigned int* h_SPD_dm(){
			if(memory_allocated && !memory_cleanup){
					return h_SPD_candidate_list_DM;
			}
			return NULL;
		}

		unsigned int* h_SPD_width(){
			if(memory_allocated && !memory_cleanup){
					return h_SPD_candidate_list_BW;
			}
			return NULL;
		}

		size_t get_SPD_nCandidates(){
			return SPD_nCandidates;
		}

		int get_current_range(){
			return (current_range>0 ? current_range-1:nRanges-1);
		}

		int get_current_tchunk(){
			return (current_range>0 ? current_time_chunk:current_time_chunk-1);
		}

		long int get_current_inc(){
			return inc;
		}

		/** \brief De-allocate memory for this pipeline instance. */
		bool cleanup() {
			if (memory_allocated && !memory_cleanup) {
				LOG(log_level::debug, "Generic Pipeline -> Memory cleanup at the end of the pipeline");

				if(do_copy_DDTR_data_to_host) {
					const int *ndms = m_ddtr_strategy.ndms_data();
					for (int i = 0; i < nRanges; i++) {
						for (int j = 0; j < ndms[i]; j++) {
							free(m_output_buffer[i][j]);
						}
						free(m_output_buffer[i]);
					}
					free(m_output_buffer);
				}

				memory_cleanup = true;
			}
			return true;
		}
	
		
	}; // class end
	

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_GENERIC_HPP

