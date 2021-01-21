#include "aa_zero_dm_outliers.hpp"

namespace astroaccelerate {

/** \brief Performs zero_dm with outliers. This is legacy code for other pipelines.*/
void zero_dm_outliers(unsigned short *const d_input, const int nchans, const int nsamp, const int nbits) {
	int threads_for_sum = 128;
	int num_blocks_t    = nsamp;
	dim3 block_size(threads_for_sum, 1);
	dim3 grid_size(num_blocks_t,1);
	float outlier_sigma = 3.0;
	cudaStream_t stream = NULL;
	int shared_memory_required = 0;
	
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
	
	
	call_kernel_zero_dm_outliers_kernel_channels(
		grid_size, 
		block_size, 
		shared_memory_required, 
		stream, 
		d_input, 
		nchans, 
		outlier_sigma, 
		nbits,
		d_normalization_factor
	);
	
	e =cudaFree(d_normalization_factor);
	if (e != cudaSuccess) {
		LOG(log_level::error, "Cannot free d_normalization_factor memory: (" + std::string(cudaGetErrorString(e)) + ")");
	}
}

/** \brief Performs zero_dm with outliers. */
void zero_dm_outliers(unsigned short *const d_input, const int nchans, const int nsamp, const int nbits, float *d_normalization_factor) {
	int threads_for_sum = 128;
	int num_blocks_t    = nsamp;
	dim3 block_size(threads_for_sum, 1);
	dim3 grid_size(num_blocks_t,1);
	float outlier_sigma = 3.0;
	cudaStream_t stream = NULL;
	int shared_memory_required = 0;
	
	call_kernel_zero_dm_outliers_kernel_channels(
		grid_size, 
		block_size, 
		shared_memory_required, 
		stream, 
		d_input, 
		nchans, 
		outlier_sigma, 
		nbits,
		d_normalization_factor
	); 
}

/** \brief Performs zero_dm with outlier rejection. */
void zero_dm_outliers_time_channels(unsigned short *const d_input, const int nchans, const int nsamp) {
	int divisions_in_t  = 224;
	int num_blocks_t    = nsamp/divisions_in_t;

	printf("\nZDM OUTLIERS!");
	printf("\n%d %d", nsamp, nchans);
	printf("\n%d %d", divisions_in_t, 1);
	printf("\n%d %d", num_blocks_t, 1);

	dim3 threads_per_block(divisions_in_t, 1);
	dim3 num_blocks(num_blocks_t,1);

	clock_t start_t, end_t;
	start_t = clock();

	call_kernel_zero_dm_outliers_kernel_one(num_blocks, threads_per_block, d_input, nchans, nsamp);
	cudaDeviceSynchronize();

	int divisions_in_c  = 100;
	int num_blocks_c    = nchans/divisions_in_c;

	printf("\nZDM OUTLIERS!");
	printf("\n%d %d", nsamp, nchans);
	printf("\n%d %d", divisions_in_c, 1);
	printf("\n%d %d", num_blocks_c, 1);

	dim3 threads_per_block_c(divisions_in_c, 1);
	dim3 c_blocks(num_blocks_c,1);
	call_kernel_zero_dm_outliers_kernel_two(c_blocks, threads_per_block_c, d_input, nchans, nsamp);
	end_t = clock();
	double time = (double)(end_t-start_t) / CLOCKS_PER_SEC;
	printf("\nPerformed ZDM: %lf (GPU estimate)", time);
}

} //namespace astroaccelerate
