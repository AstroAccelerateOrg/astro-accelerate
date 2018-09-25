#ifndef ASTRO_ACCELERATE_HOST_STRATAGY_HPP
#define ASTRO_ACCELERATE_HOST_STRATAGY_HPP

#include "host_stratagy.hpp"
#include "params.hpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "device_DDTR_plan.hpp"

class DDTR_strategy {
public:
  DDTR_Strategy(const device_DDTR_plan ddtr_plan);
  ~DDTR_Strategy();
  void setup();
  bool valid() const;
  const device_DDTR_plan& plan() const;

  //Public variables, metadata.
  //TODO: Make these variables private.
  int m_maxshift;
  int m_num_tchunks;
  int m_max_ndms;
  int m_total_ndms;
  float m_max_dm;

  // Unknowns stuff
  float m_power;

  // Sizes
  size_t m_gpu_inputsize;
  size_t m_gpu_outputsize;
  size_t m_host_inputsize;
  size_t m_host_outputsize;


  // Description of input data
  size_t m_nchans; //number of frequency channels
  size_t m_nsamp; // number of timesamples
  int m_nbits; // bit precision of input data
  float m_tsamp; //sampling time        

  
protected:
  void strategy(int *maxshift, int *max_samps, int *num_tchunks,
		int *max_ndms, int *total_ndms, float *max_dm,
		float power, int nchans, int nsamp,
		float fch1, float foff, float tsamp,
		int range, float *user_dm_low,
		float *user_dm_high, float *user_dm_step, float **dm_low,
		float **dm_high, float **dm_step, int **ndms,
		float **dmshifts, int *inBin, int ***t_processed,
		size_t *gpu_memory, int enable_analysis,
		const device_DDTR_plan &ddtr_plan);

  device_DDTR_plan m_plan;
  bool             m_ready;
};

#endif
