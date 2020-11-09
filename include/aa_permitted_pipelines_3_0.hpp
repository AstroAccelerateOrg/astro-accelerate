#ifndef ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_3_0_HPP
#define ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_3_0_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>
#include "aa_pipeline.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_ddtr_plan.hpp"


#include "aa_periodicity_strategy.hpp"

#include "aa_filterbank_metadata.hpp"
#include "aa_device_load_data.hpp"
#include "aa_bin_gpu.hpp"
#include "aa_zero_dm.hpp"
#include "aa_zero_dm_outliers.hpp"
#include "aa_corner_turn.hpp"
#include "aa_device_rfi.hpp"
#include "aa_dedisperse.hpp"
#include "aa_device_save_data.hpp"

#include "aa_gpu_timer.hpp"

#include "aa_device_periods.hpp"
#include "aa_pipeline_runner.hpp"

#include "aa_gpu_timer.hpp"

namespace astroaccelerate {

  /**
   * \class aa_permitted_pipelines_3_0 aa_permitted_pipelines_3_0.hpp "include/aa_permitted_pipelines_3_0.hpp"
   * \brief Templated class to run dedispersion and periodicity.
   * \details The class is templated over the zero_dm_type (aa_pipeline::component_option::zero_dm or aa_pipeline::component_option::zero_dm_with_outliers).
   * \author Cees Carels.
   * \date 28 November 2018.
   */

  template<aa_pipeline::component_option zero_dm_type, bool enable_old_rfi>
  class aa_permitted_pipelines_3_0 : public aa_pipeline_runner {
  public:
    aa_permitted_pipelines_3_0(const aa_ddtr_strategy &ddtr_strategy,
			       const aa_periodicity_strategy &periodicity_strategy,
			       unsigned short const*const input_buffer) {
      
    }
    
    ~aa_permitted_pipelines_3_0() {
      //Only call cleanup if memory had been allocated during setup,
      //and if the memory was not already cleaned up usingthe cleanup method.
      if(memory_allocated && !memory_cleanup) {
	cleanup();
      }      
    }

    aa_permitted_pipelines_3_0(const aa_permitted_pipelines_3_0 &) = delete;

    /** \brief Method to setup and allocate memory for the pipeline containers. */
    bool setup() override {
      if(!memory_allocated) {
	return set_data();
      }

      if(memory_allocated) {
	return true;
      }
      
      return false;
    }

    /** \brief Override base class next() method to process next time chunk. */
    bool next() override {
      if(memory_allocated) {
	aa_pipeline_runner::status tmp;
	return run_pipeline(tmp);
      }

      return false;
    }

    /** \brief Override base class next() method to process next time chunk. Also provides a status code. */
    bool next(aa_pipeline_runner::status &status_code) override {
      if(memory_allocated) {
	return run_pipeline(status_code);
      }
      
      return false;
    }
    
    /** \brief De-allocate memory for this pipeline instance. */
    bool cleanup() {
      if(memory_allocated && !memory_cleanup) {
	cudaFree(d_input);
	cudaFree(d_output);

	size_t t_processed_size = m_ddtr_strategy.t_processed().size();
	for(size_t i = 0; i < t_processed_size; i++) {
	  free(t_processed[i]);
	}
	free(t_processed);

	memory_cleanup = true;
      }
      return true;
    }
  private:
    float              ***m_output_buffer;
    int                **t_processed;
    aa_ddtr_strategy        m_ddtr_strategy;
    aa_periodicity_strategy m_periodicity_strategy;
    unsigned short     const*const m_input_buffer;
    int                num_tchunks;
    std::vector<float> dm_shifts;
    float              *dmshifts;
    int                maxshift;
    int                max_ndms;
    int                nchans;
    int                nbits;
    int                enable_zero_dm;
    int                enable_zero_dm_with_outliers;
    int                failsafe;
    long int           inc;
    float              tsamp;
    float              tsamp_original;
    int                maxshift_original;
    size_t             range;
    float              tstart_local;

    unsigned short     *d_input;
    float              *d_output;
    
    std::vector<float> dm_low;
    std::vector<float> dm_high;
    std::vector<float> dm_step;
    std::vector<int>   inBin;

    bool memory_allocated;
    bool memory_cleanup;
    bool periodicity_did_run;
    bool did_notify_of_finishing_component;
    
    //Loop counter variables
    int t;
    aa_gpu_timer       m_timer;
    
    /** \brief Allocate the GPU memory needed for dedispersion. */
    void allocate_memory_gpu(const int &maxshift, const int &max_ndms, const int &nchans, int **const t_processed, unsigned short **const d_input, float **const d_output) {

      int time_samps = t_processed[0][0] + maxshift;
      printf("\n\n\n%d\n\n\n", time_samps);
      size_t gpu_inputsize = (size_t) time_samps * (size_t) nchans * sizeof(unsigned short);

      cudaError_t e = cudaMalloc((void **) d_input, gpu_inputsize);
      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not allocate_memory_gpu cudaMalloc in aa_permitted_pipelines_3_0.hpp (" + std::string(cudaGetErrorString(e)) + ")");
      }

      size_t gpu_outputsize = 0;
      if (nchans < max_ndms) {
	gpu_outputsize = (size_t)time_samps * (size_t)max_ndms * sizeof(float);
      }
      else {
	gpu_outputsize = (size_t)time_samps * (size_t)nchans * sizeof(float);
      }

      e = cudaMalloc((void **) d_output, gpu_outputsize);
      
      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not allocate_memory_gpu cudaMalloc in aa_permitted_pipelines_3_0.hpp (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      cudaMemset(*d_output, 0, gpu_outputsize);
    }
    
    /**
     * \brief Allocate a 3D array that is an output buffer that stores dedispersed array data.
     * \details This array is used by periodicity.
     */
    void allocate_memory_cpu_output() {
      size_t estimate_outputbuffer_size = 0;
      size_t outputsize = 0;
      const size_t range = m_ddtr_strategy.get_nRanges();
      const int *ndms = m_ddtr_strategy.ndms_data();
      
      for(size_t i = 0; i < range; i++) {
        for(int j = 0; j < m_ddtr_strategy.num_tchunks(); j++) {
	  estimate_outputbuffer_size += (size_t)(t_processed[i][j]*sizeof(float)*ndms[i]);
        }
      }
      
      outputsize = 0;
      m_output_buffer = (float ***) malloc(range * sizeof(float **));
      for(size_t i = 0; i < range; i++) {
        int total_samps = 0;
        for(int k = 0; k < num_tchunks; k++) {
	  total_samps += t_processed[i][k];
	}
        m_output_buffer[i] = (float **) malloc(ndms[i] * sizeof(float *));
        for (int j = 0; j < ndms[i]; j++) {
	  m_output_buffer[i][j] = (float *) malloc(( total_samps ) * sizeof(float));
        }
        outputsize += ( total_samps ) * ndms[i] * sizeof(float);
      }
    }
    
    /** \brief Method that allocates all memory for this pipeline. */
    bool set_data() {
      num_tchunks = m_ddtr_strategy.num_tchunks();
      size_t t_processed_size = m_ddtr_strategy.t_processed().size();

      t_processed = new int*[t_processed_size];
      for(size_t i = 0; i < t_processed_size; i++) {
	t_processed[i] = new int[m_ddtr_strategy.t_processed().at(i).size()];
      }

      for(size_t i = 0; i < t_processed_size; i++) {
	for(size_t j = 0; j < m_ddtr_strategy.t_processed().at(i).size(); j++) {
	  t_processed[i][j] = m_ddtr_strategy.t_processed().at(i).at(j);
	}
      }

      dm_shifts                       = m_ddtr_strategy.dmshifts();
      dmshifts                        = dm_shifts.data();
      maxshift                        = m_ddtr_strategy.maxshift();
      max_ndms                        = m_ddtr_strategy.max_ndms();
      nchans                          = m_ddtr_strategy.metadata().nchans();
      nbits                           = m_ddtr_strategy.metadata().nbits();
      enable_zero_dm                  = 0;
      enable_zero_dm_with_outliers    = 0;
      failsafe                        = 0;
      inc                             = 0;
      tsamp                           = m_ddtr_strategy.metadata().tsamp();
      tsamp_original                  = tsamp;
      maxshift_original               = maxshift;
      range                           = m_ddtr_strategy.get_nRanges();
      tstart_local                    = 0.0;

      //Allocate GPU memory
      d_input                         = NULL;
      d_output                        = NULL;

      //Allocate GPU memory for dedispersion
      allocate_memory_gpu(maxshift, max_ndms, nchans, t_processed, &d_input, &d_output);
      //Allocate memory for CPU output for periodicity
      allocate_memory_cpu_output();
      
      //Put the dm low, high, step struct contents into separate arrays again.
      //This is needed so that the kernel wrapper functions don't need to be modified.
      dm_low.resize(m_ddtr_strategy.get_nRanges());
      dm_high.resize(m_ddtr_strategy.get_nRanges());
      dm_step.resize(m_ddtr_strategy.get_nRanges());
      inBin.resize(m_ddtr_strategy.get_nRanges());
      for(size_t i = 0; i < m_ddtr_strategy.get_nRanges(); i++) {
	dm_low[i]   = m_ddtr_strategy.dm(i).low;
	dm_high[i]  = m_ddtr_strategy.dm(i).high;
	dm_step[i]  = m_ddtr_strategy.dm(i).step;
	inBin[i]    = m_ddtr_strategy.dm(i).inBin;
      }
      memory_allocated = true;
      return true;
    }

    /**
     * \brief Run the pipeline by processing the next time chunk of data.
     * \returns A boolean to indicate whether further time chunks are available to process (true) or not (false).
     */
    bool run_pipeline(aa_pipeline_runner::status &status_code) {
      printf("NOTICE: Pipeline start/resume run_pipeline_3_0.\n");
      if(t >= num_tchunks) {

	//In order to return the status code before the next pipeline component starts,
	//the function first returns true and sets the status code to "finished_component".
	//The caller calls the pipeline again, and then the component starts.
	if(!did_notify_of_finishing_component) {
	  m_timer.Stop();
	  float time = m_timer.Elapsed() / 1000;
	  printf("\n\n === OVERALL DEDISPERSION THROUGHPUT INCLUDING SYNCS AND DATA TRANSFERS ===\n");
	  printf("\n(Performed Brute-Force Dedispersion: %g (GPU estimate)", time);
	  printf("\nAmount of telescope time processed: %f", tstart_local);
	  printf("\nNumber of samples processed: %ld", inc);
	  printf("\nReal-time speedup factor: %lf\n", ( tstart_local ) / time);
	  status_code = aa_pipeline_runner::status::finished_component;
	  did_notify_of_finishing_component = true;
	  return true;
	}

	// The status code is set after completion of periodicity, this is why a separate boolean is declared and returned after setting the status code.
	bool periodicity_return_value = periodicity();
	status_code = aa_pipeline_runner::status::finished;
	return periodicity_return_value;
      }
      else if(t == 0) {
	m_timer.Start();
      }
      printf("\nNOTICE: t_processed:\t%d, %d", t_processed[0][t], t);
      
      const int *ndms = m_ddtr_strategy.ndms_data();
      
      //checkCudaErrors(cudaGetLastError());
      load_data(-1, inBin.data(), d_input, &m_input_buffer[(long int) ( inc * nchans )], t_processed[0][t], maxshift, nchans, dmshifts, NULL, nbits);
      //checkCudaErrors(cudaGetLastError());

      if(zero_dm_type == aa_pipeline::component_option::zero_dm) {
	zero_dm(d_input, nchans, t_processed[0][t]+maxshift, nbits);
      }

      //checkCudaErrors(cudaGetLastError());

      if(zero_dm_type == aa_pipeline::component_option::zero_dm_with_outliers) {
	zero_dm_outliers(d_input, nchans, t_processed[0][t]+maxshift, nbits);
      }

      //checkCudaErrors(cudaGetLastError());

      corner_turn(d_input, d_output, nchans, t_processed[0][t] + maxshift);

      //checkCudaErrors(cudaGetLastError());

      if(enable_old_rfi) {
	printf("\nPerforming old GPU rfi...");
	rfi_gpu(d_input, nchans, t_processed[0][t]+maxshift);
      }

      //checkCudaErrors(cudaGetLastError());

      int oldBin = 1;
      for(size_t dm_range = 0; dm_range < range; dm_range++) {
	printf("\n\nNOTICE: %f\t%f\t%f\t%d\n", m_ddtr_strategy.dm(dm_range).low, m_ddtr_strategy.dm(dm_range).high, m_ddtr_strategy.dm(dm_range).step, m_ddtr_strategy.ndms(dm_range));
	printf("\nAmount of telescope time processed: %f\n", tstart_local);

	maxshift = maxshift_original / inBin[dm_range];

	cudaDeviceSynchronize();
	//checkCudaErrors(cudaGetLastError());

	load_data(dm_range, inBin.data(), d_input, &m_input_buffer[(long int) ( inc * nchans )], t_processed[dm_range][t], maxshift, nchans, dmshifts, NULL, nbits);

	//checkCudaErrors(cudaGetLastError());


	if (inBin[dm_range] > oldBin) {
	  bin_gpu(d_input, d_output, nchans, t_processed[dm_range - 1][t] + maxshift * inBin[dm_range]);
	  ( tsamp ) = ( tsamp ) * 2.0f;
	}

	//checkCudaErrors(cudaGetLastError());

	dedisperse(dm_range, t_processed[dm_range][t], inBin.data(), dmshifts, d_input, d_output, NULL, nchans, &tsamp, dm_low.data(), dm_step.data(), ndms, nbits, failsafe);

	//checkCudaErrors(cudaGetLastError());

    for (size_t k = 0; k < (size_t) ndms[dm_range]; k++) {
      size_t device_offset = (size_t) (k * (size_t) t_processed[dm_range][t]);
      size_t host_offset = (size_t) (inc / inBin[dm_range]);
      size_t data_size = (size_t) (sizeof(float) * (size_t) t_processed[dm_range][t]);
      save_data_offset(d_output, device_offset, m_output_buffer[dm_range][k], host_offset, data_size);
    }
	
	oldBin = inBin[dm_range];
      }
      
      inc = inc + t_processed[0][t];
      printf("\nNOTICE: INC:\t%ld\n", inc);
      tstart_local = ( tsamp_original * inc );
      tsamp = tsamp_original;
      maxshift = maxshift_original;

      ++t;
      printf("NOTICE: Pipeline ended run_pipeline_3_0 over chunk %d / %d.\n", t, num_tchunks);
      status_code = aa_pipeline_runner::status::has_more;
      return true;
    }

    bool periodicity() {
      if(periodicity_did_run) return false;
      cleanup();
      aa_gpu_timer timer;
      timer.Start();
      const int *ndms =	m_ddtr_strategy.ndms_data();
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
      printf("\nReal-time speedup factor: %f\n", ( tstart_local ) / ( time ));
      periodicity_did_run = true;

      for(size_t i = 0; i < range; i++) {
	for(int j = 0; j < ndms[i]; j++) {
	  free(m_output_buffer[i][j]);
	}
	free(m_output_buffer[i]);
      }
      free(m_output_buffer);
      
      return true;
    }
  };

  template<> inline aa_permitted_pipelines_3_0<aa_pipeline::component_option::zero_dm, false>::aa_permitted_pipelines_3_0(const aa_ddtr_strategy &ddtr_strategy,
															  const aa_periodicity_strategy &periodicity_strategy,
															  unsigned short const*const input_buffer) : m_ddtr_strategy(ddtr_strategy),
																				     m_periodicity_strategy(periodicity_strategy),
																				     m_input_buffer(input_buffer),
																				     memory_allocated(false),
																				     memory_cleanup(false),
																				     periodicity_did_run(false),
																				     did_notify_of_finishing_component(false),
																				     t(0) {
    
  }
  
  template<> inline aa_permitted_pipelines_3_0<aa_pipeline::component_option::zero_dm, true>::aa_permitted_pipelines_3_0(const aa_ddtr_strategy &ddtr_strategy,
															 const aa_periodicity_strategy &periodicity_strategy,
															 unsigned short const*const input_buffer) : m_ddtr_strategy(ddtr_strategy),
																				    m_periodicity_strategy(periodicity_strategy),
																				    m_input_buffer(input_buffer),
																				    memory_allocated(false),
																				    memory_cleanup(false),
																				    periodicity_did_run(false),
																				    did_notify_of_finishing_component(false),
																				    t(0) {
  }
  
  template<> inline aa_permitted_pipelines_3_0<aa_pipeline::component_option::zero_dm_with_outliers, false>::aa_permitted_pipelines_3_0(const aa_ddtr_strategy &ddtr_strategy,
																	const aa_periodicity_strategy &periodicity_strategy,
																	unsigned short const*const input_buffer) : m_ddtr_strategy(ddtr_strategy),
																						   m_periodicity_strategy(periodicity_strategy),
																						   m_input_buffer(input_buffer),
																						   memory_allocated(false),
																						   memory_cleanup(false),
																						   periodicity_did_run(false),
																						   did_notify_of_finishing_component(false),
																						   t(0) {
  }
  
  template<> inline aa_permitted_pipelines_3_0<aa_pipeline::component_option::zero_dm_with_outliers, true>::aa_permitted_pipelines_3_0(const aa_ddtr_strategy &ddtr_strategy,
																       const aa_periodicity_strategy &periodicity_strategy,
																       unsigned short const*const input_buffer) : m_ddtr_strategy(ddtr_strategy),
																						  m_periodicity_strategy(periodicity_strategy),
																						  m_input_buffer(input_buffer),
																						  memory_allocated(false),
																						  memory_cleanup(false),
																						  periodicity_did_run(false),
																						  did_notify_of_finishing_component(false),
																						  t(0) {    
  }
  
  template<> inline aa_permitted_pipelines_3_0<aa_pipeline::component_option::empty, true>::aa_permitted_pipelines_3_0(const aa_ddtr_strategy &ddtr_strategy,
														       const aa_periodicity_strategy &periodicity_strategy,
														       unsigned short const*const input_buffer) : m_ddtr_strategy(ddtr_strategy),
																				  m_periodicity_strategy(periodicity_strategy),
																				  m_input_buffer(input_buffer),
																				  memory_allocated(false),
																				  memory_cleanup(false),
																				  periodicity_did_run(false),
																				  did_notify_of_finishing_component(false),
																				  t(0) {    
  }

  template<> inline aa_permitted_pipelines_3_0<aa_pipeline::component_option::empty, false>::aa_permitted_pipelines_3_0(const aa_ddtr_strategy &ddtr_strategy,
															const aa_periodicity_strategy &periodicity_strategy,
															unsigned short const*const input_buffer) : m_ddtr_strategy(ddtr_strategy),
																				   m_periodicity_strategy(periodicity_strategy),
																				   m_input_buffer(input_buffer),
																				   memory_allocated(false),
																				   memory_cleanup(false),
																				   periodicity_did_run(false),
																				   did_notify_of_finishing_component(false),
																				   t(0) {
  }
  
} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_3_0_HPP
