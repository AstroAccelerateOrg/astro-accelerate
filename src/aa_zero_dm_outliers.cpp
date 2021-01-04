#include "aa_zero_dm_outliers.hpp"

namespace astroaccelerate {

/** \brief Performs zero_dm. */
void zero_dm_outliers(unsigned short *const d_input, const int nchans, const int nsamp, const int nbits) {
	int threads_for_sum = 128;
	int num_blocks_t    = nsamp;
	dim3 block_size(threads_for_sum, 1);
	dim3 grid_size(num_blocks_t,1);
	float normalization_factor = ((pow(2,nbits)-1)/2);
	float outlier_sigma = 3.0;
	cudaStream_t stream = NULL;
	int shared_memory_required = 0;
	
	//printf("------> Zero DM with outliers INFO\n");
	//printf("block_size=[%d;%d;%d]; grid_size=[%d;%d;%d]\n", block_size.x, block_size.y, block_size.z, grid_size.x, grid_size.y, grid_size.z);
	//printf("Normalization factor=%f; outlier_sigma=%f;\n", normalization_factor, outlier_sigma);
	
	call_kernel_zero_dm_outliers_kernel_channels(
		grid_size, 
		block_size, 
		shared_memory_required, 
		stream, 
		d_input, 
		nchans, 
		outlier_sigma, 
		nbits,
		normalization_factor
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
