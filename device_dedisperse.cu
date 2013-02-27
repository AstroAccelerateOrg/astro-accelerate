
#define NOPSMIN 3.6
#define NOPSMAX 4.6

extern "C" void dedisperse(size_t inputsize, float *d_input, size_t outputsize, float *d_output, int nchans, int nsamp, int maxshift, float dm_low, int ndms, int kernel_type, float tsamp, float dm_step);

//{{{ dedisperse 

void dedisperse(size_t inputsize, float *d_input, size_t outputsize, float *d_output, int nchans, int nsamp, int maxshift, float dm_low, int ndms, int kernel_type, float tsamp, float dm_step) {

	//{{{ Set the timing parameters

//	cudaEvent_t start, stop;
//	float time[10], sum;
//	int i;

//	for(i = 0; i < 10; i++) time[i] = 0.0;

//	i = 0;
	
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);

	//}}}
	
	//{{{ Dedisperse data on the GPU 

	float startdm = dm_low;

	int num_reg         = NUMREG;
	int divisions_in_t  = DIVINT;
	int divisions_in_dm = DIVINDM;
	int num_blocks_t    = (nsamp-maxshift)/(divisions_in_t * num_reg);
	int num_blocks_dm   = ndms/divisions_in_dm;

	dim3 threads_per_block(divisions_in_t, divisions_in_dm);
	dim3 num_blocks(num_blocks_t,num_blocks_dm);

//	printf("\n\tndms:\t\t%d", ndms);
//	printf("\n\tdm_step:\t%f", dm_step);
//	printf("\n\tnsamp:\t\t%d", nsamp);
//	printf("\n\ttsamp:\t\t%lf", tsamp);
//	printf("\n\tnum_acc:\t%d", num_reg);
//	printf("\n\tBlocksize(x,y):\t%d,%d",  divisions_in_t, divisions_in_dm);
//	printf("\n\tGridsize(x,y):\t%d,%d\n", num_blocks_t, num_blocks_dm);

//	printf("\n\tkernelStart"), fflush(stdout);
//	cudaEventRecord(start,0);

	if(kernel_type == 0) {
		printf("\n\tUsing Shared Memory Algorithm");
		cudaFuncSetCacheConfig(shared_dedisperse_loop, cudaFuncCachePreferShared);
		shared_dedisperse_loop<<< num_blocks, threads_per_block >>>(d_output, d_input, (float)(startdm/tsamp), (float)(dm_step/tsamp));
	} else if(kernel_type == 1) {
		printf("\n\tUsing L1 cache Algorithm");
		cudaFuncSetCacheConfig(cache_dedisperse_loop, cudaFuncCachePreferL1);
		cache_dedisperse_loop<<<  num_blocks, threads_per_block >>>(d_output, d_input, (float)(startdm/tsamp), (float)(dm_step/tsamp));
	} else if(kernel_type == 2) {
		printf("\n\tUsing contiguous L1 cache Algorithm");
		cudaFuncSetCacheConfig(cache_contiguous_loop, cudaFuncCachePreferL1);
		cache_contiguous_loop<<<  num_blocks, threads_per_block >>>(d_output, d_input, (float)(startdm/tsamp), (float)(dm_step/tsamp));
	}  else if(kernel_type == 3) {
		printf("\n\tUsing contiguous Shared Memory Algorithm");
		cudaFuncSetCacheConfig(shared_contiguous_loop, cudaFuncCachePreferShared);
		shared_contiguous_loop<<< num_blocks, threads_per_block >>>(d_output, d_input, (float)(startdm/tsamp), (float)(dm_step/tsamp));
	}

//	cudaEventRecord(stop, 0);
//	cudaEventSynchronize(stop);
//	cudaEventElapsedTime(&time[i], start, stop);
//	printf("\n\tkernelStop"), fflush(stdout);
//	printf("\n\tPerformed Brute-Force Dedispersion:\t\t%f ms (GPU estimate)", time[i]);

//	printf("\n\n\tReal-time speedup factor:\t\t\t%f", ((nsamp - maxshift)*tsamp)/(time[i]/1000));
//	printf("\n\tGops based on %.2f ops per channel per tsamp:\t%f",NOPSMIN,((1.0*NOPSMIN*ndms*nchans*(nsamp - maxshift))/(time[i]/1000))/1000000000);
//	printf("\n\tGops based on %.2f ops per channel per tsamp:\t%f",NOPSMAX,((1.0*NOPSMAX*ndms*nchans*(nsamp - maxshift))/(time[i]/1000))/1000000000);
//	i++;

	//}}}


//	for(i = 0; i < 10; i++) sum += time[i];

//	printf("\n\n\tGPU total time:\t\t\t\t\t%lf s", sum/1000);

}

//}}}

