//#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <cufft.h>
#include "AstroAccelerate/params.h"
#include "device_stats_kernel.cu"
#include "helper_cuda.h"

//{{{ Return stats 

void stats_gpu(cudaEvent_t event, cudaStream_t stream, int samps, float *mean, float *stddev, float *h_signal_power, float *d_signal_power)
{

	int a, j;
	int trials = ( 2 * ACCMAX + ACCSTEP ) / ACCSTEP;
	// int chunk = omp_get_num_procs();
	/*
	 // Calculate the mean
	 double	total=0.0;
	 #pragma omp parallel for default(shared) private(a,j) schedule(static,chunk) reduction(+:total)
	 for(a = 0; a < trials; a++) {
	 for(j=0;j< (samps/2);j++){
	 total += ((double)(h_signal_power[j+(a)*(samps/2)]));
	 //		if(h_signal_power[j+(a)*(samps/2)] > 2465757) printf("\n%d %d %f", a, j, h_signal_power[j+(a)*(samps/2)]);
	 }
	 }
	 *mean = (float)(total/(double)((samps/2)*(trials+1)));  // Mean for data sample

	 // Calculate standard deviation
	 total = 0.0;
	 #pragma omp parallel for default(shared) private(a,j) schedule(static,chunk) reduction(+:total)
	 for(a = 0; a < trials; a++) {
	 for(j=0; j<(samps/2); j++){
	 total += (double)((h_signal_power[j+(a)*(samps/2)]-(*mean))*(h_signal_power[j+(a)*(samps/2)]-(*mean)));
	 }
	 }
	 *stddev = (float)sqrt(abs(total) / (double)((samps/2)*(trials+1))); // Stddev for data sample
	 */
	int half_samps = samps / 2;
	int acc_size = half_samps * trials;

	int divisions = STATST;
	int blocks = (int) floor((float) acc_size / divisions / STATSLOOP);

	dim3 threads_per_block(divisions);
	dim3 num_blocks(blocks);

	int size = (int) floor((float) acc_size / STATSLOOP);

	float* d_sum;
	checkCudaErrors(cudaMalloc((void** )&d_sum, size * sizeof(float)));

	float* d_sum_square;
	checkCudaErrors(cudaMalloc((void** )&d_sum_square, size * sizeof(float)));

	float* h_sum;
	checkCudaErrors(cudaMallocHost((void** )&h_sum, size * sizeof(float)));

	float* h_sum_square;
	checkCudaErrors(cudaMallocHost((void** )&h_sum_square, size * sizeof(float)));

	cudaStreamWaitEvent(stream, event, 0);
	stats_kernel<<<num_blocks, threads_per_block, 0, stream>>>(half_samps, d_sum, d_sum_square, d_signal_power);
	getLastCudaError("power_kernel failed");
	cudaEventRecord(event, stream);

	cudaStreamWaitEvent(stream, event, 0);
	checkCudaErrors(cudaMemcpyAsync(h_sum, d_sum, size * sizeof(float), cudaMemcpyDeviceToHost, stream));
	checkCudaErrors(cudaMemcpyAsync(h_sum_square, d_sum_square, size * sizeof(float), cudaMemcpyDeviceToHost, stream));
	cudaEventRecord(event, stream);
	cudaStreamSynchronize(stream);

	float total_sum = 0.0;
	float total_sum_square = 0.0;
//#pragma omp parallel for default(shared) private(a) schedule(static,chunk) reduction(+:total_sum,total_sum_square)
	for (a = 0; a < size; a++)
	{
		total_sum += ( h_sum[a] );
		total_sum_square += ( h_sum_square[a] );
	}
	*mean = (float) ( total_sum / ( acc_size ) );  // Mean for data sample
	*stddev = (float) sqrt(( total_sum_square - acc_size * ( *mean ) * ( *mean ) ) / ( acc_size - 1 ));
	//printf("\nM:\t%f, S:\t%f", *mean, *stddev);

	cudaFree(d_sum);
	cudaFree(d_sum_square);
	cudaFreeHost(h_sum);
	cudaFreeHost(h_sum_square);

}

//}}}

