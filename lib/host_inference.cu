/* Written by Andrew Parker as part of COSIT.
 * COSIT was a NVIDIA funded project and supervised by W Armour
 * 07/2013-10/2013
 */ 

#include <stdio.h>

void cpu_blocked_bootstrap(float *data_array, unsigned int num_els, unsigned int num_bins, unsigned int num_boots, float *mean_data_extern, float *cpu_mean_extern, float *cpu_sd_extern){

	//double cpu_timer, timeSinceLastTimer;

	float partial_squaresum=0.0f;
	float cpu_mean=0.0f, cpu_sd=0.0f, mean_data=0.0f;
	
	//cpu_timer=omp_get_time_of_day();

	//bin array will contain the partial sum of many numbers
	float *bin_array;
	bin_array = (float*) malloc (num_bins*sizeof(float));
	for(int j=0;j<num_bins;j++) bin_array[j]=0.0f;

	//bootstrap mean will contain the bootstrap resamples
	float *bootstrap_mean;
	bootstrap_mean = (float*) malloc(num_boots*sizeof(float));
	for(int j = 0; j<num_boots; j++) bootstrap_mean[j]=0.0f;

	//care should be taken here when num_bins*num_boots is large
	int *rand_array;
	rand_array = (int*) malloc (num_bins*num_boots*sizeof(int));
	for(int j=0;j<num_bins*num_boots;j++) rand_array[j]=0;

	//timeSinceLastTimer = elapsed_time(&cpu_timer);
	//printf("CPU Time to assign memory:\t%f \n",timeSinceLastTimer);

	//rand array contains a list of discrete random numbers between 0 and binsize-1
	for(int i = 0; i<num_bins*num_boots;i++) {
        	rand_array[i]= num_bins*(rand()/(float) RAND_MAX);
	}

	//timeSinceLastTimer = elapsed_time(&cpu_timer);
	//printf("CPU Time to create random array:\t%f \n",timeSinceLastTimer);

	//bin array contains the sum of the numbers in each bin
	for(int j = 0; j<num_bins; j++){
        	for(int i = 0; i< (num_els/num_bins); i++){
			bin_array[j]+= data_array[i+j*(num_els/num_bins)];
        	}
	}

	//timeSinceLastTimer = elapsed_time(&cpu_timer);
	//printf("CPU Time to create bin array:\t%f \n",timeSinceLastTimer);


	//bootstrap array contains an array of bootstrap means
	//for each bootstrap, sum the num_bins sums into bootstrap_array; then divide these sums by num_elements
	for(int j = 0; j<num_boots; j++){
        	for(int i = 0; i< num_bins; i++){
			bootstrap_mean[j]+=bin_array[rand_array[i+j*num_bins]];
		}
        }

	for(int j=0; j<num_boots; j++) bootstrap_mean[j]/=num_els;

	//cpu_calcs = elapsed_time(&cpu_timer);
	//printf("CPU Time to compute bootstrap means:\t%f \n",cpu_calcs);

	//compute mean of bootstrap means
	for(int j = 0; j<num_boots; j++){
        	cpu_mean += bootstrap_mean[j];
	}
	cpu_mean/=num_boots;

	// compute sd of the bootstrap means
	for(int j = 0; j<num_boots; j++){
		partial_squaresum += (bootstrap_mean[j]-cpu_mean)*(bootstrap_mean[j]-cpu_mean);
	}
	cpu_sd=sqrt((partial_squaresum)/(num_boots-1));

	// compute mean of data
	for (int i = 0; i<num_bins; i++){
        	mean_data += bin_array[i];
	}
	mean_data/= num_els;

	//timeSinceLastTimer = elapsed_time(&cpu_timer);
	//printf("CPU Time to calculate statistics:\t%f \n",timeSinceLastTimer);
	//printf("\nCPU data mean = %f\n",mean_data);
	//printf("CPU bootstrap_mean = %f\n",cpu_mean);
	//printf("CPU bootstrap_sd = %f\n",cpu_sd);
	//printf("\n%f\t%f\t%f\t%f", dm, mean_data, cpu_mean, cpu_sd);
	*mean_data_extern=mean_data;
	*cpu_mean_extern=cpu_mean;
	*cpu_sd_extern=cpu_sd;

	//free data
        free(bin_array);
        free(bootstrap_mean);
        free(rand_array);
}
