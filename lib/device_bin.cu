//#include <helper_cuda.h>
// #include <omp.h>
#include <time.h>
#include <stdio.h>
#include "AstroAccelerate/params.h"
#include "device_binning_kernel.cu"

//extern "C" void bin_gpu(float *bin_buffer, float *input_buffer, int nchans, int nsamp);

void bin_gpu(unsigned short *d_input, float *d_output, int nchans, int nsamp)
{

	int divisions_in_t = BINDIVINT;
	int divisions_in_f = BINDIVINF;
	int num_blocks_t = (int) ( ( nsamp + 1 ) / ( 2 * divisions_in_t ) );
	int num_blocks_f = nchans / divisions_in_f;

	dim3 threads_per_block(divisions_in_t, divisions_in_f);
	dim3 num_blocks(num_blocks_t, num_blocks_f);
	//printf("\nDIVISOR:\t%f", (float)(nsamp)/(2*divisions_in_t));
	//printf("\ndivisions_in_t:%d\tdivisions_in_f:%d",divisions_in_t, divisions_in_f);
	//printf("\nnum_blocks_t:%d\tnum_blocks_f:%d\n",num_blocks_t,num_blocks_f);

	//cudaFuncSetCacheConfig(bin, cudaFuncCachePreferL1);
	// clock_t start_t, end_t;
	// start_t = clock();

	bin<<<num_blocks, threads_per_block>>>(d_input, d_output, nsamp);
	//	getLastCudaError("Kernel execution failed");

	int swap_divisions_in_t = CT;
	int swap_divisions_in_f = CF;
	int swap_num_blocks_t = nsamp / swap_divisions_in_t;
	int swap_num_blocks_f = nchans / swap_divisions_in_f;

	dim3 swap_threads_per_block(swap_divisions_in_t, swap_divisions_in_f);
	dim3 swap_num_blocks(swap_num_blocks_t, swap_num_blocks_f);

	cudaDeviceSynchronize();
	swap<<<swap_num_blocks, swap_threads_per_block>>>(d_input, d_output, nchans, nsamp);
	cudaDeviceSynchronize();

	// end_t = clock();
	// double time = double ( end_t - start_t ) / CLOCKS_PER_SEC;
	//printf("\nPerformed Bin: %f (GPU estimate)", time);
	//printf("\nGops based on %.2f ops per channel per tsamp: %f",14.0,((15.0*(divisions_in_t*divisions_in_f*num_blocks_t*num_blocks_f))/(time))/1000000000.0);
	//printf("\nBN Device memory bandwidth in GB/s: %f", (2*(sizeof(float)+sizeof(unsigned short))*(divisions_in_t*divisions_in_f*num_blocks_t*num_blocks_f))/(time)/1000000000.0);

	cudaMemset(d_output, 0, nchans * nsamp * sizeof(float));
	//cudaMemcpy(input_buffer, bin_buffer, sizeof(float)*nchans*(nsamp/2), cudaMemcpyDeviceToDevice);
	//	getLastCudaError("Kernel execution failed");
}

