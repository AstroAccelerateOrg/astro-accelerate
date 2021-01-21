#include "aa_zero_dm.hpp"

namespace astroaccelerate {

  /** \brief Performs zero_dm. This is legacy code for other pipelines. */
  void zero_dm(unsigned short *const d_input, const int nchans, const int nsamp, const int nbits) {
    
    int threads_for_sum  = 256;
    int num_blocks_t    = nsamp;
    int shared_memory = threads_for_sum*4;
    
    dim3 threads_per_block(threads_for_sum, 1);
    dim3 num_blocks(num_blocks_t, 1);
    
    float normalization_factor = ((pow(2,nbits)-1)/2);
	std::vector<float> local_bandpass_normalization;
	local_bandpass_normalization.resize(nchans,normalization_factor);
	float *d_normalization_factor = NULL;
	cudaError_t e;
	e = cudaMalloc((void **)d_normalization_factor, nchans*sizeof(float));
	if (e != cudaSuccess) {
		LOG(log_level::error, "Could not allocate memory for d_normalization_factor (" + std::string(cudaGetErrorString(e)) + ")");
	}
	cudaMemcpy(d_normalization_factor, local_bandpass_normalization.data(), nchans*sizeof(float), cudaMemcpyHostToDevice);
	
    call_kernel_zero_dm_kernel(num_blocks, threads_per_block, shared_memory, d_input, nchans, nsamp, nbits, d_normalization_factor);
    cudaDeviceSynchronize();
	
	e =cudaFree(d_normalization_factor);
	if (e != cudaSuccess) {
		LOG(log_level::error, "Cannot free d_normalization_factor memory: (" + std::string(cudaGetErrorString(e)) + ")");
	}
  }

  /** \brief Performs zero_dm. */
  void zero_dm(unsigned short *const d_input, const int nchans, const int nsamp, const int nbits, float *d_normalization_factor) {
    
    int threads_for_sum  = 256;
    int num_blocks_t    = nsamp;
    int shared_memory = threads_for_sum*4;
    
    dim3 threads_per_block(threads_for_sum, 1);
    dim3 num_blocks(num_blocks_t, 1);
    
    call_kernel_zero_dm_kernel(num_blocks, threads_per_block, shared_memory, d_input, nchans, nsamp, nbits, d_normalization_factor);
    cudaDeviceSynchronize();
  }
} //namespace astroaccelerate
