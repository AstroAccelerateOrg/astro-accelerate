#ifndef ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_GENERIC_HPP
#define ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_GENERIC_HPP

#define PIPELINE_ERROR_NO_ERROR 0
#define PIPELINE_ERROR_DDTR_GPU_MEMORY_FAIL 1
#define PIPELINE_ERROR_SPDT_GPU_MEMORY_FAIL 2
#define PIPELINE_ERROR_HOST_MEMORY_FAIL 3
#define PIPELINE_ERROR_GPU_FREE_MEMORY_FAIL 4
#define PIPELINE_ERROR_GENERAL_GPU_ERROR 5
#define PIPELINE_ERROR_ZERO_DM 6
#define PIPELINE_ERROR_CORNER_TURN 7
#define PIPELINE_ERROR_RFI 8
#define PIPELINE_ERROR_BINNING 9
#define PIPELINE_ERROR_DEDISPERSION 10
#define PIPELINE_ERROR_COPY_TO_HOST 11
#define PIPELINE_ERROR_COPY_TO_DEVICE 12
#define PIPELINE_ERROR_SPDT_ERROR 13


#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "aa_pipeline.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_ddtr_plan.hpp"
#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"
#include "aa_periodicity_strategy.hpp"
#include "aa_fdas_strategy.hpp"
#include "aa_jerk_strategy.hpp"

#include "aa_filterbank_metadata.hpp"
#include "aa_device_load_data.hpp"
#include "aa_bin_gpu.hpp"
#include "aa_zero_dm.hpp"
#include "aa_zero_dm_outliers.hpp"
#include "aa_corner_turn.hpp"
#include "aa_device_rfi.hpp"
#include "aa_dedisperse.hpp"
#include "aa_device_save_data.hpp"

#include "aa_device_analysis.hpp"
#include "aa_device_periods.hpp"
#include "aa_device_acceleration_fdas.hpp"
#include "aa_device_jerk_search.hpp"
#include "aa_pipeline_runner.hpp"

#include "aa_gpu_timer.hpp"
#include "aa_timelog.hpp"
#include "aa_host_cpu_msd.hpp"

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
		aa_jerk_strategy                     m_jerk_strategy;
		
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
		int                pipeline_error;

		unsigned short     *d_DDTR_input;
		float              *d_DDTR_output;
		float              *d_dm_shifts;
		float              *d_bandpass_normalization;

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
		bool do_jerk;
		
		bool do_copy_DDTR_data_to_host;
		
		bool memory_allocated;
		bool memory_cleanup;
		bool periodicity_did_run;
		bool acceleration_did_run;
		bool jerk_did_run;
		bool did_notify_of_finishing_component;

		//Loop counter variables
		int current_time_chunk;
		int current_range;
		TimeLog time_log;
		aa_gpu_timer m_timer;
		aa_gpu_timer m_local_timer;
		aa_gpu_timer m_ddtr_total_timer;
		
		float  *m_d_MSD_workarea = NULL;
		float  *m_d_MSD_interpolated = NULL;
		ushort *m_d_SPDT_output_taps = NULL;

		bool cleanup_DDTR() {
			LOG(log_level::debug, "DDTR -> Memory cleanup after de-dispersion");
			
			cudaError_t e;
			if(d_DDTR_input!=NULL) {
				e = cudaFree(d_DDTR_input);
				if (e != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_GPU_FREE_MEMORY_FAIL;
					LOG(log_level::error, "Cannot free memory (" + std::string(cudaGetErrorString(e)) + ")");
				}
				else {
					d_DDTR_input = NULL;
				}
			}
			if(d_DDTR_output!=NULL) {
				e = cudaFree(d_DDTR_output);
				if (e != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_GPU_FREE_MEMORY_FAIL;
					LOG(log_level::error, "Cannot free memory (" + std::string(cudaGetErrorString(e)) + ")");
				}
				else {
					d_DDTR_output = NULL;
				}
			}
			if(nchans>8192 && nbits!=4 && d_dm_shifts!=NULL) {
				e =cudaFree(d_dm_shifts);
				if (e != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_GPU_FREE_MEMORY_FAIL;
					LOG(log_level::error, "Cannot free memory (" + std::string(cudaGetErrorString(e)) + ")");
				}
				else {
					d_dm_shifts = NULL;
				}
			}
			if(nchans>4096 && nbits==4 && d_dm_shifts!=NULL) {
				e =cudaFree(d_dm_shifts);
				if (e != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_GPU_FREE_MEMORY_FAIL;
					LOG(log_level::error, "Cannot free memory (" + std::string(cudaGetErrorString(e)) + ")");
				}
				else {
					d_dm_shifts = NULL;
				}
			}
			
			e = cudaFree(d_bandpass_normalization);
			if (e != cudaSuccess) {
				pipeline_error = PIPELINE_ERROR_GPU_FREE_MEMORY_FAIL;
				LOG(log_level::error, "Cannot free d_bandpass_normalization memory: (" + std::string(cudaGetErrorString(e)) + ")");
			}
			else {
				d_bandpass_normalization = NULL;
			}
            
			//----------------- Single pulse memory de-allocation --->
			if(do_single_pulse_detection){
				if(m_d_MSD_workarea!=NULL) {
					e = cudaFree(m_d_MSD_workarea);
					if (e != cudaSuccess) {
						pipeline_error = PIPELINE_ERROR_GPU_FREE_MEMORY_FAIL;
						LOG(log_level::error, "Cannot free memory (" + std::string(cudaGetErrorString(e)) + ")");
					}
					else {
						m_d_MSD_workarea = NULL;
					}
				}
				if(m_d_SPDT_output_taps!=NULL) {
					e = cudaFree(m_d_SPDT_output_taps);
					if (e != cudaSuccess) {
						pipeline_error = PIPELINE_ERROR_GPU_FREE_MEMORY_FAIL;
						LOG(log_level::error, "Cannot free memory (" + std::string(cudaGetErrorString(e)) + ")");
					}
					else {
						m_d_SPDT_output_taps = NULL;
					}
				}
				if(m_d_MSD_interpolated!=NULL) {
					e = cudaFree(m_d_MSD_interpolated);
					if (e != cudaSuccess) {
						pipeline_error = PIPELINE_ERROR_GPU_FREE_MEMORY_FAIL;
						LOG(log_level::error, "Cannot free memory (" + std::string(cudaGetErrorString(e)) + ")");
					}
					else {
						m_d_MSD_interpolated = NULL;
					}
				}

					
				if(h_SPD_candidate_list_DM!=NULL)  {
					free(h_SPD_candidate_list_DM);
					h_SPD_candidate_list_DM = NULL;
				}
				if(h_SPD_candidate_list_TS!=NULL) {
					free(h_SPD_candidate_list_TS);
					h_SPD_candidate_list_TS = NULL;
				}
				if(h_SPD_candidate_list_SNR!=NULL) {
					free(h_SPD_candidate_list_SNR);
					h_SPD_candidate_list_SNR = NULL;
				}
				if(h_SPD_candidate_list_BW!=NULL) {
					free(h_SPD_candidate_list_BW);
					h_SPD_candidate_list_BW = NULL;
				}
					
			}
				
			/*
			if (memory_allocated && !memory_cleanup) {
				LOG(log_level::debug, "DDTR -> Memory cleanup after de-dispersion");
				cudaFree(d_DDTR_input);
				cudaFree(d_DDTR_output);
				if( (nchans > 8192) && (nbits != 4) ) cudaFree(d_dm_shifts);
				if( (nbits == 4) && (nchans > 4096) ) cudaFree(d_dm_shifts);
				
				if(do_single_pulse_detection){
					cudaFree(m_d_MSD_workarea);
					cudaFree(m_d_SPDT_output_taps);
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
			*/
			
			if(pipeline_error!=0){
				return false;
			}
			else return true;
		}
		
		/** \brief Allocate the GPU memory needed for dedispersion. */
		bool allocate_gpu_memory_DDTR(){
				//const int &maxshift, const int &max_ndms, const int &nchans, int **const t_processed, unsigned short **const d_DDTR_input, float **const d_DDTR_output, float **const d_dm_shifts) {
			
			size_t free_memory = 0, total_memory = 0, required_memory = 0;
			cudaMemGetInfo(&free_memory,&total_memory);
			int time_samps = t_processed[0][0] + maxshift;
			size_t gpu_inputsize = (size_t)time_samps * (size_t)nchans * sizeof(unsigned short);
			size_t gpu_outputsize = 0;
			if (nchans < max_ndms) {
				gpu_outputsize = (size_t)time_samps * (size_t)max_ndms * sizeof(float);
			}
			else {
				gpu_outputsize = (size_t)time_samps * (size_t)nchans * sizeof(float);
			}
			
			//---------------------------------> DEBUG info
			printf("\n\nMemory allocation for the DDTR\n");
			printf("Number of time samples: %d;\n", t_processed[0][0]);
			printf("maxshift: %d\n", maxshift);
			printf("Number of time samples + maxshift: %d;\n", time_samps);
			printf("Number of channels: %d\n", nchans);
			printf("Maximum number of DM-trials: %d;\n", max_ndms);
			printf("Available memory: %zu bytes = %0.3f MB;\n", free_memory, ((double) free_memory)/(1024.0*1024.0));
			printf("DDTR input size: %zu bytes = %0.3f MB;\n", gpu_inputsize, ((double) gpu_inputsize)/(1024.0*1024.0));
			printf("DDTR output size:  %zu bytes = %0.3f MB;\n", gpu_outputsize, ((double) gpu_outputsize)/(1024.0*1024.0));
			if(nchans>8192) {
				printf("Channel shifts: %zu bytes = %0.3f MB;\n", nchans*sizeof(float), ((double) nchans*sizeof(float))/(1024.0*1024.0));
				required_memory = required_memory + nchans*sizeof(float);
			}
			required_memory = required_memory + gpu_outputsize + gpu_inputsize;
			printf("Total memory required: %zu bytes = %0.3f MB;\n", required_memory, ((double) required_memory)/(1024.0*1024.0));
			
			//---------------------------------> Allocations
			cudaError_t e;
			e = cudaMalloc((void **)&d_DDTR_input, gpu_inputsize);
			if (e != cudaSuccess) {
				pipeline_error = PIPELINE_ERROR_DDTR_GPU_MEMORY_FAIL;
				LOG(log_level::error, "Could not allocate memory for d_DDTR_input using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}

			e = cudaMalloc((void **)&d_DDTR_output, gpu_outputsize);
			if (e != cudaSuccess) {
				pipeline_error = PIPELINE_ERROR_DDTR_GPU_MEMORY_FAIL;
				LOG(log_level::error, "Could not allocate memory for d_DDTR_output using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}
			cudaMemset(d_DDTR_output, 0, gpu_outputsize);

			if( (nchans>8192) && (nbits != 4) ){		
				e = cudaMalloc((void **) &d_dm_shifts, nchans*sizeof(float));
				if (e != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_DDTR_GPU_MEMORY_FAIL;
					LOG(log_level::error, "Could not allocate memory for d_dm_shifts using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
				}
			}

			if( (nchans>4096) && (nbits==4) ){
				e = cudaMalloc((void **) &d_dm_shifts, nchans*sizeof(float));
				if (e != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_DDTR_GPU_MEMORY_FAIL;
					LOG(log_level::error, "Could not allocate memory for d_dm_shifts (4-bit) using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
				}
			}
			
			//Allocation of bandpass normalization for zerodm
			e = cudaMalloc((void **) &d_bandpass_normalization, nchans*sizeof(float));
			if (e != cudaSuccess) {
				pipeline_error = PIPELINE_ERROR_DDTR_GPU_MEMORY_FAIL;
				LOG(log_level::error, "Could not allocate memory for d_bandpass_normalization using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}
			
			if(pipeline_error!=PIPELINE_ERROR_NO_ERROR) return false;
			else return true;
		}


		/** \brief Allocate memory for single pulse detection (SPD) on the host. */
		void allocate_cpu_memory_SPD(size_t max_nDMs, size_t timesamples_per_chunk){
			SPD_max_peak_size        = (size_t)(max_nDMs*timesamples_per_chunk/2);
			
			h_SPD_candidate_list_DM  = (unsigned int*)malloc(SPD_max_peak_size*sizeof(unsigned int));
			h_SPD_candidate_list_TS  = (unsigned int*)malloc(SPD_max_peak_size*sizeof(unsigned int));
			h_SPD_candidate_list_SNR = (float*)malloc(SPD_max_peak_size*sizeof(float));
			h_SPD_candidate_list_BW  = (unsigned int*)malloc(SPD_max_peak_size*sizeof(unsigned int));
			if(h_SPD_candidate_list_DM==NULL || h_SPD_candidate_list_TS==NULL || h_SPD_candidate_list_SNR==NULL || h_SPD_candidate_list_BW==NULL) {
				pipeline_error = PIPELINE_ERROR_HOST_MEMORY_FAIL;
				LOG(log_level::error, "Could not allocate memory on the host for single pulse detection candidates");
			}
		}
		
		
		/** \brief Allocate memory for single pulse detection (SPD) on the device.	*/
		void allocate_gpu_memory_SPD(float **const d_MSD_workarea, unsigned short **const d_SPDT_output_taps, float **const d_MSD_interpolated, const unsigned long int &MSD_maxtimesamples, const size_t &MSD_profile_size) {
			cudaError_t e;
			
			e = cudaMalloc((void **)d_MSD_workarea, MSD_maxtimesamples*5.5*sizeof(float));
			if (e != cudaSuccess) {
				pipeline_error = PIPELINE_ERROR_SPDT_GPU_MEMORY_FAIL;
				LOG(log_level::error, "Could not allocate memory for d_MSD_workarea using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}

			e = cudaMalloc((void **) &(*d_SPDT_output_taps), sizeof(ushort)*2*MSD_maxtimesamples);
			if (e != cudaSuccess) {
				pipeline_error = PIPELINE_ERROR_SPDT_GPU_MEMORY_FAIL;
				LOG(log_level::error, "Could not allocate memory for d_SPDT_output_taps using cudaMalloc in aa_permitted_pipelines_generic.hpp (" + std::string(cudaGetErrorString(e)) + ")");
			}

			e = cudaMalloc((void **)d_MSD_interpolated, sizeof(float)*MSD_profile_size);
			if (e != cudaSuccess) {
				pipeline_error = PIPELINE_ERROR_SPDT_GPU_MEMORY_FAIL;
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

			if(do_copy_DDTR_data_to_host){
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
			}
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
			d_DDTR_input  = NULL;
			d_DDTR_output = NULL;
			d_dm_shifts   = NULL;
			
			//Data for zero_dm
			d_bandpass_normalization = NULL;
			
			//Allocate GPU memory for dedispersion
			//allocate_gpu_memory_DDTR(maxshift, max_ndms, nchans, t_processed, &d_DDTR_input, &d_DDTR_output, &d_dm_shifts);
			allocate_gpu_memory_DDTR();
			
			if(m_ddtr_strategy.bandpass_normalization_size() == (size_t) nchans){
				cudaMemcpy(d_bandpass_normalization, m_ddtr_strategy.bandpass_normalization_pointer(), nchans*sizeof(float), cudaMemcpyHostToDevice);
			}
			else pipeline_error=PIPELINE_ERROR_ZERO_DM;
			
			
			//Allocate GPU memory for SPD (i.e. analysis)
			if(do_single_pulse_detection && pipeline_error==0) {
				allocate_cpu_memory_SPD(max_ndms, t_processed[0][0]);
				allocate_gpu_memory_SPD(&m_d_MSD_workarea, &m_d_SPDT_output_taps, &m_d_MSD_interpolated, m_analysis_strategy.MSD_data_info(), m_analysis_strategy.MSD_profile_size_in_bytes());
			}
			//Allocate memory for CPU output for periodicity
			if(pipeline_error==0) {
				allocate_memory_cpu_output();
			}
			
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
			
			if(pipeline_error!=0) {
				
				memory_allocated = false;
				return false;
			}
			else {
				memory_allocated = true;
				return true;
			}
		}

		bool input_data_preprocessing(){
			cudaError_t cuda_error;
			const aa_pipeline::component_option opt_set_bandpass_average = aa_pipeline::component_option::set_bandpass_average;
			if (m_pipeline_options.find(opt_set_bandpass_average) != m_pipeline_options.end()) {
				LOG(log_level::debug, "Calculating zero DM bandpass (average)... ");
				m_local_timer.Start();
				
				float *h_new_bandpass;
				h_new_bandpass = new float[nchans];
				size_t nsamples = m_ddtr_strategy.metadata().nsamples();

				
				// Place code here
				call_cpu_msd(h_new_bandpass, m_input_buffer, (size_t)nchans, nsamples);
				
				cuda_error = cudaMemcpy(d_bandpass_normalization, h_new_bandpass, nchans*sizeof(float), cudaMemcpyHostToDevice);
				if (cuda_error != cudaSuccess) {
					pipeline_error=PIPELINE_ERROR_ZERO_DM;
				}
				delete[] h_new_bandpass;
				
				m_local_timer.Stop();
				time_log.adding("PREP", "Bandpass", m_local_timer.Elapsed());
			}
			
			if(pipeline_error!=0) {
				return false;
			}
			else {
				return true;
			}
		}

		/** \brief Method that pipeline flags according to requested pipeline components. */
		void set_pipeline_flags(){
			//----> Components
			//const aa_pipeline::component cmp_dedispersion = aa_pipeline::component::dedispersion;
			const aa_pipeline::component cmp_analysis     = aa_pipeline::component::analysis;
			const aa_pipeline::component cmp_periodicity  = aa_pipeline::component::periodicity;
			const aa_pipeline::component cmp_fdas         = aa_pipeline::component::fdas;
			const aa_pipeline::component cmp_jerk         = aa_pipeline::component::jerk;
			
			//----> Component options
			const aa_pipeline::component_option opt_copy_ddtr_data_to_host = aa_pipeline::component_option::copy_ddtr_data_to_host;
			
			do_dedispersion = true; // default component
			if(m_pipeline_components.find(cmp_analysis) != m_pipeline_components.end()) do_single_pulse_detection = true; 
			else do_single_pulse_detection = false;
			if(m_pipeline_components.find(cmp_periodicity) != m_pipeline_components.end()) do_periodicity_search = true;
			else do_periodicity_search = false;
			if(m_pipeline_components.find(cmp_fdas) != m_pipeline_components.end()) do_fdas = true;
			else do_fdas = false;
			if(m_pipeline_components.find(cmp_jerk) != m_pipeline_components.end()) do_jerk = true;
			else do_jerk = false;
			
			
			do_copy_DDTR_data_to_host = false;
			if(m_pipeline_options.find(opt_copy_ddtr_data_to_host) != m_pipeline_options.end()) do_copy_DDTR_data_to_host = true;
			if(do_periodicity_search) do_copy_DDTR_data_to_host = true;
			if(do_fdas) do_copy_DDTR_data_to_host = true;
			if(do_jerk) do_copy_DDTR_data_to_host = true;
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
			const aa_pipeline::component_option opt_input_DDTR_normalization  = aa_pipeline::component_option::input_DDTR_normalization;
			const aa_pipeline::component_option opt_output_DDTR_normalization = aa_pipeline::component_option::output_DDTR_normalization;
			const aa_pipeline::component_option opt_old_rfi                = aa_pipeline::component_option::old_rfi;
			
			//------------------------------------> Error checking
			if(pipeline_error!=0) {
				// TODO: add error strings
				LOG(log_level::error, "Pipeline encountered an error.");
			}
			
			cudaError_t CUDA_error;
			CUDA_error = cudaGetLastError();
			if(CUDA_error != cudaSuccess) {
				pipeline_error = PIPELINE_ERROR_GENERAL_GPU_ERROR;
				LOG(log_level::error, "GPU error at the pipeline start. (" + std::string(cudaGetErrorString(CUDA_error)) + ")");
			}
			//-------------------------------------<
			
			
			printf("NOTICE: Pipeline start/resume run_pipeline_generic.\n");
			if (current_time_chunk >= num_tchunks) {
				//------------------> End of DDTR + SPD
				if (!did_notify_of_finishing_component) {
					m_timer.Stop();
					time_log.adding("Total", "total", m_timer.Elapsed());
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
				
				
				//------------------> JERK search
				if (!jerk_did_run && do_jerk) {
					bool jerk_return_value = jerk_search();
					status_code = aa_pipeline_runner::status::finished;
					return jerk_return_value;
				}
				//--------------------------------------------------------------------------------<

				return false; // In this case, there are no more chunks to process, and periodicity and acceleration both ran.
			}
			else if (current_time_chunk == 0 && current_range == 0) {
				m_timer.Start();
			}
			
			//----------------------------------------------------------------------->
			//------------------> Dedispersion
			const int *ndms = m_ddtr_strategy.ndms_data();
			
			if(current_range==0 && pipeline_error==PIPELINE_ERROR_NO_ERROR) {
				m_ddtr_total_timer.Start();
				printf("\nNOTICE: t_processed:\t%d, %d", t_processed[0][current_time_chunk], current_time_chunk);

				//---------> Load chunks
				m_local_timer.Start();
				load_chunk_data(d_DDTR_input, &m_input_buffer[(long int)(inc * nchans)], t_processed[0][current_time_chunk], maxshift_original, nchans, dmshifts, d_dm_shifts, nbits);
				m_local_timer.Stop();
				time_log.adding("DDTR", "Host_To_Device", m_local_timer.Elapsed());
				
				CUDA_error = cudaGetLastError();
				if(CUDA_error != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_COPY_TO_DEVICE;
					LOG(log_level::error, "GPU error at ZeroDM kernel. (" + std::string(cudaGetErrorString(CUDA_error)) + ")");
				}
				//-------------------------<
				
				
				//---------> Zero DM
				if (m_pipeline_options.find(opt_zero_dm) != m_pipeline_options.end()) {
					LOG(log_level::debug, "Performing zero DM...");
					m_local_timer.Start();
					zero_dm(d_DDTR_input, nchans, t_processed[0][current_time_chunk] + maxshift_original, nbits, d_bandpass_normalization);
					m_local_timer.Stop();
					time_log.adding("DDTR", "Zero_DM", m_local_timer.Elapsed());
				}
				else if (m_pipeline_options.find(opt_zero_dm_with_outliers) != m_pipeline_options.end()) {
					LOG(log_level::debug, "Performing zero DM with outliers...");
					m_local_timer.Start();
					zero_dm_outliers(d_DDTR_input, nchans, t_processed[0][current_time_chunk] + maxshift_original, nbits, d_bandpass_normalization);
					m_local_timer.Stop();
					time_log.adding("DDTR", "Zero_DM_outliers", m_local_timer.Elapsed());
				}
				
				CUDA_error = cudaGetLastError();
				if(CUDA_error != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_ZERO_DM;
					LOG(log_level::error, "GPU error at ZeroDM kernel. (" + std::string(cudaGetErrorString(CUDA_error)) + ")");
				}
				//-------------------------<

	
				//---------> Input DDTR data normalization
				if (m_pipeline_options.find(opt_input_DDTR_normalization) != m_pipeline_options.end()) {
					LOG(log_level::debug, "Performing input DDTR normalization...");
					m_local_timer.Start();
					size_t nTimesamples = t_processed[0][current_time_chunk] + maxshift_original;
					zero_dm_normalization_dm(d_DDTR_input, nchans, nTimesamples, nbits);
					m_local_timer.Stop();
					time_log.adding("DDTR", "input_DDTR_norm", m_local_timer.Elapsed());
				}
				//-------------------------<
				//---------> Corner turn
				m_local_timer.Start();
				corner_turn(d_DDTR_input, d_DDTR_output, nchans, t_processed[0][current_time_chunk] + maxshift_original);
				m_local_timer.Stop();
				time_log.adding("DDTR", "Corner_Turn", m_local_timer.Elapsed());
				
				CUDA_error = cudaGetLastError();
				if(CUDA_error != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_CORNER_TURN;
					LOG(log_level::error, "GPU error at Corner turn. (" + std::string(cudaGetErrorString(CUDA_error)) + ")");
				}
				//-------------------------<


				//---------> old RFI
				if (m_pipeline_options.find(opt_old_rfi) != m_pipeline_options.end()) {
					printf("\nPerforming old GPU rfi...");
					m_local_timer.Start();
					rfi_gpu(d_DDTR_input, nchans, t_processed[0][current_time_chunk] + maxshift_original);
					m_local_timer.Stop();
					time_log.adding("DDTR", "RFI_GPU", m_local_timer.Elapsed());
					
					CUDA_error = cudaGetLastError();
					if(CUDA_error != cudaSuccess) {
						pipeline_error = PIPELINE_ERROR_RFI;
						LOG(log_level::error, "GPU error at RFI. (" + std::string(cudaGetErrorString(CUDA_error)) + ")");
					}
				}
				//-------------------------<

				oldBin = 1;
				m_ddtr_total_timer.Stop();
				time_log.adding("DDTR", "total", m_ddtr_total_timer.Elapsed());
			}

			
			if(current_time_chunk<num_tchunks && pipeline_error==PIPELINE_ERROR_NO_ERROR){
				m_ddtr_total_timer.Start();
				int dm_range = current_range;
				float tsamp = tsamp_original*((float) inBin[dm_range]);
				printf("\n\nNOTICE: %f\t%f\t%f\t%d\n", m_ddtr_strategy.dm(dm_range).low, m_ddtr_strategy.dm(dm_range).high, m_ddtr_strategy.dm(dm_range).step, m_ddtr_strategy.ndms(dm_range));
				printf("\nAmount of telescope time processed: %f\n", tstart_local);
				
				maxshift = maxshift_original / inBin[dm_range];

				cudaDeviceSynchronize();
				m_local_timer.Start();
				set_dedispersion_constants(t_processed[dm_range][current_time_chunk], maxshift);
				m_local_timer.Stop();
				time_log.adding("DDTR", "Host_To_Device",m_local_timer.Elapsed());
				
				CUDA_error = cudaGetLastError();
				if(CUDA_error != cudaSuccess) {
					pipeline_error = PIPELINE_ERROR_COPY_TO_HOST;
					LOG(log_level::error, "GPU error at Dedispersion. (" + std::string(cudaGetErrorString(CUDA_error)) + ")");
				}

				
				if (inBin[dm_range] > oldBin) {
					m_local_timer.Start();
					bin_gpu(d_DDTR_input, d_DDTR_output, nchans, t_processed[dm_range - 1][current_time_chunk] + maxshift * inBin[dm_range]);
					m_local_timer.Stop();
					time_log.adding("DDTR", "Binning",m_local_timer.Elapsed());
					
					CUDA_error = cudaGetLastError();
					if(CUDA_error != cudaSuccess) {
						pipeline_error = PIPELINE_ERROR_BINNING;
						LOG(log_level::error, "GPU error at Binning. (" + std::string(cudaGetErrorString(CUDA_error)) + ")");
					}
				}
				
				int kernel_error;				
				m_local_timer.Start();
				kernel_error = dedisperse(dm_range, t_processed[dm_range][current_time_chunk], inBin.data(), dmshifts, d_DDTR_input, d_DDTR_output, d_dm_shifts, nchans, &tsamp, dm_low.data(), dm_step.data(), ndms, nbits, failsafe);
				m_local_timer.Stop();
				time_log.adding("DDTR","Dedispersion",m_local_timer.Elapsed());

				if(kernel_error != 0) {
					pipeline_error = PIPELINE_ERROR_DEDISPERSION;
					//LOG(log_level::error, "GPU error at Dedispersion. (" + std::string(cudaGetErrorString(CUDA_error)) + ")");
					LOG(log_level::error, "GPU error at Dedispersion.");
				}
				
				//---------> Output DDTR data normalization
				if (m_pipeline_options.find(opt_output_DDTR_normalization) != m_pipeline_options.end()) {
					LOG(log_level::debug, "Performing output DDTR normalization...");
					m_local_timer.Start();
					size_t nTimesamples = t_processed[dm_range][current_time_chunk];
					size_t nDMs = ndms[dm_range];
					post_DDTR_normalization(d_DDTR_output, nTimesamples, nDMs);
					m_local_timer.Stop();
					time_log.adding("DDTR", "output_DDTR_norm", m_local_timer.Elapsed());
				}
				//-------------------------<
				
				//-----------> Copy data to the host
				if(do_copy_DDTR_data_to_host){
					m_local_timer.Start();
						save_data_offset_stream(dm_range, current_time_chunk, t_processed, inc, inBin.data(), ndms, d_DDTR_output, m_output_buffer);
					m_local_timer.Stop();
					time_log.adding("DDTR", "Device_To_Host", m_local_timer.Elapsed());
				}
				m_ddtr_total_timer.Stop();
				time_log.adding("DDTR", "total", m_ddtr_total_timer.Elapsed());
				//--------------------------------------------<

				

				//------------------> Single pulse detection
				if (do_single_pulse_detection && pipeline_error==PIPELINE_ERROR_NO_ERROR) {
					const bool dump_to_disk = true;
					const bool dump_to_user = false;
					bool SPDT_no_error;
					analysis_output output;
					SPD_nCandidates = 0;
					SPDT_no_error = analysis_GPU(
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
						m_d_SPDT_output_taps,
						m_d_MSD_interpolated,
						m_analysis_strategy.MSD_data_info(),
						m_analysis_strategy.enable_msd_baseline_noise(),
						dump_to_disk,
						dump_to_user,
						output);
					printf("Number of candidates: %zu\n", SPD_nCandidates);
					if(SPDT_no_error==false) {
						pipeline_error = PIPELINE_ERROR_SPDT_ERROR;
					}
				}
				//--------------------------------------------------------------------------------<
				oldBin = inBin[dm_range];
			}

			printf("NOTICE: Pipeline ended run_pipeline_generic over chunk %d / %d and range %d / %d.\n", current_time_chunk, num_tchunks, current_range, nRanges);
			
			++current_range;
			if(current_range>=nRanges) {
				inc = inc + t_processed[0][current_time_chunk];
				tstart_local = (tsamp_original * inc);
				printf("\nNOTICE: INC:\t%ld\n", inc);
				++current_time_chunk;
				current_range = 0;
			}
			
			
			// Anomalous pipeline termination
			if(pipeline_error!=PIPELINE_ERROR_NO_ERROR) {
				cleanup_DDTR();
				status_code = aa_pipeline_runner::status::error;
				return false;
			}
			status_code = aa_pipeline_runner::status::has_more;
			return true;
		} // run_pipeline

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
			time_log.adding("Periodicity", "total", timer.Elapsed());
			time_log.adding("Total", "total", timer.Elapsed());

			//float time = timer.Elapsed()/1000;
			// printf("\n\n === OVERALL PERIODICITY THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");
			// printf("\nPerformed Periodicity Location: %f (GPU estimate)", time);
			// printf("\nAmount of telescope time processed: %f", tstart_local);
			// printf("\nNumber of samples processed: %ld", inc);
			// printf("\nReal-time speedup factor: %f\n", (tstart_local) / (time) );
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
			time_log.adding("FDAS", "total", timer.Elapsed());
			time_log.adding("Total", "total", timer.Elapsed());
			float time = timer.Elapsed()/1000;
			printf("\n\n === OVERALL TDAS THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");

			printf("\nPerformed Acceleration Location: %lf (GPU estimate)", time);
			printf("\nAmount of telescope time processed: %f", tstart_local);
			printf("\nNumber of samples processed: %ld", inc);
			printf("\nReal-time speedup factor for FDAS: %lf\n", (tstart_local) / (time));
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
		
		bool jerk_search() {
			const int *ndms = m_ddtr_strategy.ndms_data();
			int nRanges = m_ddtr_strategy.get_nRanges();
			
			aa_gpu_timer timer;
			timer.Start();
			
			jerk_search_from_ddtr_plan(m_output_buffer, m_jerk_strategy, dm_low.data(), dm_step.data(), ndms, tsamp_original, inBin.data(), nRanges);
			
			timer.Stop();
			time_log.adding("JERK", "total", timer.Elapsed());
			time_log.adding("Total", "total", timer.Elapsed());
			
			jerk_did_run = true;
			return (true);
		}

	public:
		aa_permitted_pipelines_generic(
		const aa_pipeline::pipeline &pipeline_components,
		const aa_pipeline::pipeline_option &pipeline_options,
		const aa_ddtr_strategy &ddtr_strategy,
		const aa_analysis_strategy &analysis_strategy,
		const aa_periodicity_strategy &periodicity_strategy,
		const aa_fdas_strategy &fdas_strategy,
		const aa_jerk_strategy &jerk_strategy,
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
		m_jerk_strategy(jerk_strategy),
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
		jerk_did_run(false),
		did_notify_of_finishing_component(false),
		current_time_chunk(0),
		current_range(0),
		m_d_MSD_workarea(NULL),
		m_d_MSD_interpolated(NULL),
		m_d_SPDT_output_taps(NULL) {
	}

		~aa_permitted_pipelines_generic() {
			//Only call cleanup if memory had been allocated during setup,
			//and if the memory was not already cleaned up usingthe cleanup method.
			if (memory_allocated && !memory_cleanup) {
				cleanup();
			}
			if(pipeline_error==PIPELINE_ERROR_NO_ERROR){
				time_log.print(tstart_local);
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
			pipeline_error = PIPELINE_ERROR_NO_ERROR;
			set_pipeline_flags();
			
			if (!memory_allocated) {
				set_data();
				if(memory_allocated) input_data_preprocessing();
			}

			if (memory_allocated) {
				return true;
			}
			else {
				cleanup_DDTR();
				return false;
			}
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
			int temp_range = (current_range>0 ? current_range-1:nRanges-1);
			if(temp_range==nRanges-1 && current_time_chunk>0){
				return(inc-t_processed[0][current_time_chunk-1]);
			}
			else {
				return(inc);
			}
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
			
			if (m_pipeline_options.find(aa_pipeline::component_option::timelog_export_to_file) != m_pipeline_options.end()) {
				time_log.print_to_file();
			}

			return true;
		}
		
	}; // class end

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_GENERIC_HPP


