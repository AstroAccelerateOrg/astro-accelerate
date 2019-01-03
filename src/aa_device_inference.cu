/* Written by Andrew Parker as part of COSIT.
 * COSIT was a NVIDIA funded project and supervised by W Armour
 * 07/2013-10/2013
 */

#include <stdio.h>
#include <cuda.h>
#include <cufft.h>
#include <curand.h>

namespace astroaccelerate {

#define BSTRAP 2
#define ELPB 32
#define NTDM 4

  //__constant__ int num_els;
  //__constant__ float sum;
  //__constant__ float mean;
  //__constant__ int bins;

  //void cpu_blocked_bootstrap(float *data_array, unsigned int num_els, unsigned int num_bins, unsigned int num_boots, float *mean_data_extern, float *cpu_mean_extern, float *cpu_sd_extern){

  // ************************************* GPU Blocked Bootstrap ************************************** //

  // ************************************* Bootstrap Kernel ************************************** //
  // most bias, but quickest (rand in (0,bins-blockDim.x))
  __global__ void bootstrap(int bins, int num_els, int num_boots, float *g_idata, double *g_odata, unsigned int *g_irand)
  {
    float myResample = 0.0f;

    unsigned int constant = ( 4294967295 / ( bins - blockDim.x ) );
    int constant2 = blockIdx.x * bins;
    int dmid = bins * ( blockDim.y * blockIdx.y + threadIdx.y );
    for (int i = 0; i < bins; i++)
      {
	int rid = ( g_irand[constant2 + i] / ( constant ) );
	myResample += g_idata[rid + dmid + threadIdx.x];
      }
    dmid = num_boots * ( blockDim.y * blockIdx.y + threadIdx.y );
    g_odata[dmid + threadIdx.x + blockDim.x * blockIdx.x] = ( (double) myResample / (double) num_els );
  }

  //some bias, slower than bootstrap (rand in (0,bins-1) but with modulo to loop)
  __global__ void bootstrap2(int bins, int num_els, int num_boots, float *g_idata, double *g_odata, unsigned int *g_irand)
  {
    float myResample = 0.0f;

    int constant = ( 4294967295 / ( bins ) );
    int constant2 = blockIdx.x * bins;
    int dmid = bins * ( blockDim.y * blockIdx.y + threadIdx.y );
    for (int i = 0; i < bins; i++)
      {

	int rid = g_irand[constant2 + i] / constant;

	myResample += g_idata[dmid + ( ( rid + threadIdx.x ) % bins )];
      }
    dmid = num_boots * ( blockDim.y * blockIdx.y + threadIdx.y );
    g_odata[dmid + threadIdx.x + blockDim.x * blockIdx.x] = ( (double) myResample / (double) num_els );
  }

  //least amount of bias, slower than bootstrap2 (rand array is much bigger, each thread chooses a sample independent of other threads).
  __global__ void bootstrap3(int bins, int num_els, int num_boots, float *g_idata, double *g_odata, unsigned int *g_irand)
  {
    float myResample;

    int constant = ( 4294967295 / ( bins ) );
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int dmid = bins * ( blockDim.y * blockIdx.y + threadIdx.y );
    for (int i = 0; i < bins; i++)
      {

	int rid = g_irand[id * bins + i] / constant;

	myResample += g_idata[dmid + rid];
      }
    dmid = num_boots * ( blockDim.y * blockIdx.y + threadIdx.y );
    g_odata[dmid + threadIdx.x + blockDim.x * blockIdx.x] = ( (double) myResample / (double) num_els );

  }

  void gpu_blocked_bootstrap(float **d_idata, int dms_to_average, int num_els, int ndms, int num_bins, int num_boots, double *mean_boot_out, double *mean_data_out, double *sd_boot_out)
  {

    int local_ndms = ndms / dms_to_average;

    float *d_bin_array;
    cudaMalloc((void**) &d_bin_array, local_ndms * num_bins * sizeof(float));
    float *bin_array;
    bin_array = (float *) malloc(local_ndms * num_bins * sizeof(float));
    float *random_bin_array;
    random_bin_array = (float *) malloc(local_ndms * num_bins * sizeof(float));

    // Do I need to worry about rounding error??
    int nsamp_per_bin = num_els / num_bins;
    for (int dm_count = 0; dm_count < local_ndms; dm_count++)
      {
	for (int b = 0; b < num_bins; b++)
	  {
	    float local_sum = 0.0f;
	    for (int a = 0; a < dms_to_average; a++)
	      {
		for (int j = b * nsamp_per_bin; j < ( b + 1 ) * nsamp_per_bin; j++)
		  {
		    local_sum += d_idata[dm_count * dms_to_average + a][j];
		  }
	      }
	    bin_array[dm_count * num_bins + b] = local_sum / dms_to_average;
	  }
      }

    /*
    //Our Array exclusive_random is used to assigns subbin sums into random bins
    int **exclusive_random;
    exclusive_random=(int **)malloc(local_ndms*sizeof(int *));
    for(int dm_count=0; dm_count<local_ndms; dm_count++) {
    exclusive_random[dm_count]=(int *) malloc(num_bins*sizeof(int));
    for(int i = 0; i < num_bins; i++) exclusive_random[dm_count][i]=i;
    }

    for(int dm_count=0; dm_count<local_ndms; dm_count++) {
    //Fischer shuffle algorithm to replace our bins randomly
    //sample without replacement once for indexes to replace
    for(int i=0;i<num_bins-1;i++){
    int Rand=(((num_bins-i)*rand())/(float) RAND_MAX );
    int temp=exclusive_random[dm_count][num_bins-i-1];
    exclusive_random[dm_count][num_bins-i-1]=exclusive_random[dm_count][Rand];
    exclusive_random[dm_count][Rand]=temp;
    }
    }

    for(int dm_count=0; dm_count<local_ndms; dm_count++) {
    for (int j=0; j<num_bins; j++) {
    random_bin_array[dm_count*num_bins + exclusive_random[dm_count][j]]=bin_array[dm_count*num_bins +j];
    }
    }
    cudaMemcpy(d_bin_array, random_bin_array, local_ndms*num_bins*sizeof(float), cudaMemcpyHostToDevice);

    for(int dm_count=0; dm_count<local_ndms; dm_count++) {
    free(exclusive_random[dm_count]);
    }
    free(exclusive_random);
    */
    //Our Array exclusive_random is used to assigns subbin sums into random bins
    int *exclusive_random;
    exclusive_random = (int *) malloc(num_bins * sizeof(int));
    for (int i = 0; i < num_bins; i++)
      exclusive_random[i] = i;

    //Fischer shuffle algorithm to replace our bins randomly
    //sample without replacement once for indexes to replace
    for (int i = 0; i < num_bins - 1; i++)
      {
	int Rand = ( ( ( num_bins - i ) * rand() ) / (float) RAND_MAX );
	int temp = exclusive_random[num_bins - i - 1];
	exclusive_random[num_bins - i - 1] = exclusive_random[Rand];
	exclusive_random[Rand] = temp;
      }

    for (int dm_count = 0; dm_count < local_ndms; dm_count++)
      {
	for (int j = 0; j < num_bins; j++)
	  {
	    random_bin_array[dm_count * num_bins + exclusive_random[j]] = bin_array[dm_count * num_bins + j];
	  }
      }
    cudaMemcpy(d_bin_array, random_bin_array, local_ndms * num_bins * sizeof(float), cudaMemcpyHostToDevice);
    free(exclusive_random);

    // create array of random numbers: we need num_subboots*num_bins random numbers
    unsigned int *d_irand;
    curandGenerator_t gen;
    int num_subboots = num_boots / ELPB;
    if (BSTRAP == 3)
      {
	( cudaMalloc((void**) &d_irand, ( num_bins ) * ( num_boots ) * sizeof(unsigned int)) );
      }
    else
      {
	( cudaMalloc((void**) &d_irand, ( num_bins ) * ( num_subboots ) * sizeof(unsigned int)) );
      }
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);

    if (BSTRAP == 3)
      {
	curandGenerate(gen, d_irand, ( num_bins ) * ( num_boots ));
      }
    else
      {
	curandGenerate(gen, d_irand, ( num_bins ) * ( num_subboots ));
      }
    cudaDeviceSynchronize();

    int divisions_in_dm = NTDM;
    int num_blocks_dm = local_ndms / divisions_in_dm;
    //printf("\n%f",(float)local_ndms/divisions_in_dm);

    int divisions_in_b;
    int num_blocks_b;

    divisions_in_b = ELPB;
    num_blocks_b = num_subboots;

    dim3 threads_per_block_boot(divisions_in_b, divisions_in_dm);
    dim3 num_blocks_boot(num_blocks_b, num_blocks_dm);

    double *d_odata;
    cudaMalloc((void**) &d_odata, local_ndms * num_boots * sizeof(double));

    if (BSTRAP == 1)
      bootstrap<<<num_blocks_boot, threads_per_block_boot>>>(num_bins, num_els, num_boots, d_bin_array, d_odata, d_irand);
    if (BSTRAP == 2)
      bootstrap2<<<num_blocks_boot, threads_per_block_boot>>>(num_bins, num_els, num_boots, d_bin_array, d_odata, d_irand);
    if (BSTRAP == 3)
      bootstrap3<<<num_blocks_boot, threads_per_block_boot>>>(num_bins, num_els, num_boots, d_bin_array, d_odata, d_irand);
    //cudaCheckMsg("bootstrap kernel execution failed");

    double *boots_array;
    boots_array = (double*) malloc(local_ndms * ( num_boots ) * sizeof(double));
    for (int i = 0; i < num_boots * local_ndms; i++)
      boots_array[i] = 0.0;

    cudaDeviceSynchronize();
    cudaMemcpy(boots_array, d_odata, local_ndms * num_boots * sizeof(double), cudaMemcpyDeviceToHost);

    //finally calculate the mean and standard deviation for the bootstraps
    for (int dm_count = 0; dm_count < local_ndms; dm_count++)
      {
	double sd_boots = 0.0f, mean_data = 0.0f, mean_boots = 0.0f;
	for (int i = 0; i < num_boots; i++)
	  {
	    mean_boots += boots_array[dm_count * num_boots + i];
	  }
	mean_boots /= (double) num_boots;
	mean_boot_out[dm_count] = mean_boots;

	for (int i = 0; i < num_bins; i++)
	  {
	    mean_data += bin_array[dm_count * num_bins + i];
	  }
	mean_data /= (double) num_els;
	mean_data_out[dm_count] = mean_data;

	for (int i = 0; i < num_boots; i++)
	  {
	    sd_boots += ( boots_array[dm_count * num_boots + i] - mean_boots ) * ( boots_array[dm_count * num_boots + i] - mean_boots );
	  }
	sd_boots = sqrt(( sd_boots / (double) ( num_boots - 1 ) ));
	sd_boot_out[dm_count] = sd_boots;
      }

    ( cudaFree(d_bin_array) );
    ( cudaFree(d_odata) );
    ( cudaFree(d_irand) );

    free(bin_array);
    free(random_bin_array);
    free(boots_array);
  }

} //namespace astroaccelerate
