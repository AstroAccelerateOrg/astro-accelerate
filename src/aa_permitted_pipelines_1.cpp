//
//  aa_permitted_pipelines_1.cpp
//  aapipeline
//
//  Created by Cees Carels on Friday 02/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_permitted_pipelines_1.hpp"

namespace astroaccelerate {

inline void save_data_offset(float *device_pointer, int device_offset, float *host_pointer, int host_offset, size_t size) {
    cudaMemcpy(host_pointer + host_offset, device_pointer + device_offset, size, cudaMemcpyDeviceToHost);
}
  
  bool aa_permitted_pipelines_1::run_pipeline(std::vector<float> &output_buffer) {
  printf("NOTICE: Pipeline start/resume run_pipeline_1.\n");
  const int *ndms = m_ddtr_strategy.ndms_data();
  bool dump_ddtr_output = true;
  if(t >= num_tchunks) return false;//In this case, there are no more chunks to process.
  
  printf("\nNOTICE: t_processed:\t%d, %d", t_processed[0][t], t);
  
  checkCudaErrors(cudaGetLastError());
  load_data(-1, inBin.data(), d_input, &m_input_buffer[(long int) ( inc * nchans )], t_processed[0][t], maxshift, nchans, dmshifts);
  checkCudaErrors(cudaGetLastError());
  
  if (enable_zero_dm) {
    zero_dm(d_input, nchans, t_processed[0][t]+maxshift, nbits);
  }
  
  checkCudaErrors(cudaGetLastError());
  
  
  if (enable_zero_dm_with_outliers) {
    zero_dm_outliers(d_input, nchans, t_processed[0][t]+maxshift);
  }
  
  checkCudaErrors(cudaGetLastError());
  
  corner_turn(d_input, d_output, nchans, t_processed[0][t] + maxshift);
  
  checkCudaErrors(cudaGetLastError());
  
  int oldBin = 1;
  for(size_t dm_range = 0; dm_range < range; dm_range++) {
    printf("\n\nNOTICE: %f\t%f\t%f\t%d\n", m_ddtr_strategy.dm(dm_range).low, m_ddtr_strategy.dm(dm_range).high, m_ddtr_strategy.dm(dm_range).step, m_ddtr_strategy.ndms(dm_range));
    printf("\nAmount of telescope time processed: %f\n", tstart_local);
    
    maxshift = maxshift_original / inBin[dm_range];
    
    cudaDeviceSynchronize();
    checkCudaErrors(cudaGetLastError());
    
    load_data(dm_range, inBin.data(), d_input, &m_input_buffer[(long int) ( inc * nchans )], t_processed[dm_range][t], maxshift, nchans, dmshifts);
    
    checkCudaErrors(cudaGetLastError());
    
    
    if (inBin[dm_range] > oldBin) {
      bin_gpu(d_input, d_output, nchans, t_processed[dm_range - 1][t] + maxshift * inBin[dm_range]);
      ( tsamp ) = ( tsamp ) * 2.0f;
    }
    
    checkCudaErrors(cudaGetLastError());
    
    dedisperse(dm_range, t_processed[dm_range][t], inBin.data(), dmshifts, d_input, d_output, nchans, &tsamp, dm_low.data(), dm_step.data(), ndms, nbits, failsafe);

    if(dump_ddtr_output) {
      //Resize vector to contain the output array
      size_t total_samps = 0;
      for (int k = 0; k < num_tchunks; k++) {
	total_samps += t_processed[dm_range][k];
      }
      output_buffer.resize(total_samps);
      for (int k = 0; k < ndms[dm_range]; k++) {
	save_data_offset(d_output, k * t_processed[dm_range][t], output_buffer.data(), inc / inBin[dm_range], sizeof(float) * t_processed[dm_range][t]);
      }
    }
    
    checkCudaErrors(cudaGetLastError());
    
    oldBin = inBin[dm_range];
  }
  
  inc = inc + t_processed[0][t];
  printf("\nNOTICE: INC:\t%ld\n", inc);
  tstart_local = ( tsamp_original * inc );
  tsamp = tsamp_original;
  maxshift = maxshift_original;

  ++t;
  printf("NOTICE: Pipeline ended run_pipeline_1 over chunk %d / %d.\n", t, num_tchunks);
  return true;
}
  
} //namespace astroaccelerate
