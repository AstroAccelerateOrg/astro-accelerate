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
  DDTR_Strategy(const device_DDTR_plan& ddtr_plan);
  ~DDTR_Strategy();
  void setup();
  bool valid() const;
  const device_DDTR_plan& plan() const;
  int totalNumberOfTimeChunks() const;
  int totalNumberOfTimeChunks() const;

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

  float *m_dmshifts; //size of nchans
  size_t **t_processed;

protected:
  int Allocate_dmshifts();
  int Allocate_t_processed_outer();

  void strategy(int *max_samps, float fch1, float foff, float tsamp,
		int *max_ndms, int *total_ndms, float *max_dm,
		int nchans, int nsamp,
		size_t *gpu_memory, int enable_analysis);

  device_DDTR_plan m_plan;
  bool             m_ready;
};

#endif
