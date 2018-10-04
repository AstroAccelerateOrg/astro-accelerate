#ifndef ASTRO_ACCELERATE_HOST_STRATAGY_HPP
#define ASTRO_ACCELERATE_HOST_STRATAGY_HPP

#include "host_stratagy.hpp"
#include "InputDataMeta.hpp"
#include "params.hpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "device_DDTR_plan.hpp"

class DDTR_strategy {
  public:
    DDTR_Strategy(device_DDTR_plan const& ddtr_plan, InputDataMeta const& input_data, size_t gpu_memory);
    ~DDTR_Strategy();
    void setup();
    bool valid() const;

    const device_DDTR_plan& plan() const;
    const InputDataMeta& data_meta() const;

    int totalNumberOfTimeChunks() const;
    int totalNumberOfTimeChunks() const;

    //Public variables, metadata.
    //TODO: Make these variables private.
    int m_maxshift;
    int m_num_tchunks;
    int m_max_ndms;
    int m_total_ndms;
    float m_max_dm;
    std::size_t m_gpu_memory;

    // Unknowns stuff
    float m_power;

    // Sizes
    size_t m_gpu_inputsize;
    size_t m_gpu_outputsize;
    size_t m_host_inputsize;
    size_t m_host_outputsize;


    float *m_dmshifts; //size of nchans
    size_t **t_processed;
    int* m_ndms;
    int m_total_ndms;

    // Description of input data
    InputDataMeta m_meta;

  protected:
    // protected members
    int strategy(int *max_ndms, int *total_ndms, int enable_analysis);

  private:
    int Allocate_dmshifts();
    int Allocate_t_processed_outer();

  private:
    device_DDTR_plan m_plan;
    bool             m_ready;
};

#endif
