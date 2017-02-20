//Added by Karel Adamek

#include "AstroAccelerate/params.h"
#include "device_SPS_long_kernel.cu"


void Get_next_iteration_parameters(int old_nBoxcars, int new_nBoxcars, int nDMs, int nTimesamples, int *nBlocks, int *unprocessed_samples, int *decimated_timesamples, int *shift, int *output_shift, int *startTaps, int  *iteration) {
	int t_decimated_timesamples;
	int t_nBlocks, t_nRest, t_Elements_per_block;

	*shift = old_nBoxcars/2;
	if(old_nBoxcars==0) *output_shift = 0;
	else *output_shift = *output_shift + nDMs*(nTimesamples>>(*iteration));
	*startTaps = *startTaps + old_nBoxcars*(1<<(*iteration));
	
	if(old_nBoxcars==0) *iteration = 0;
	else *iteration = *iteration + 1;
	
	t_decimated_timesamples = nTimesamples>>(*iteration);
	t_Elements_per_block=PD_NTHREADS*2-new_nBoxcars;
	t_nBlocks=(t_decimated_timesamples - new_nBoxcars + 1)/t_Elements_per_block;
	t_nRest=(t_decimated_timesamples - new_nBoxcars + 1) - t_nBlocks*t_Elements_per_block;
	
	*decimated_timesamples = t_decimated_timesamples;
	*unprocessed_samples = t_nRest + new_nBoxcars-1;
	*nBlocks = t_nBlocks;
}

int Get_max_iteration(int max_boxcar_width, int *PD_plan, int PD_plan_size){
	int startTaps, iteration;
	
	startTaps = 0;
	for(int f=0; f<PD_plan_size; f++){
		startTaps = startTaps + PD_plan[f]*(1<<f);
		if(startTaps>=max_boxcar_width) {
			iteration = f+1;
			break;
		}
	}
	
	return(iteration);
}

void PD_SEARCH_LONG_init() {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeEightByte);
}

int PD_SEARCH_LONG(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int max_boxcar_width, int nDMs, int nTimesamples, int *t_max_iterarion) {
	//---------> Task specific
	int nBlocks, nRest, Elements_per_block, unprocessed_samples;
	int shift, output_shift, iteration, max_iteration, startTaps, decimated_timesamples;
	int nBoxcars;
	int PD_plan[10]={16,16,16,16,8,8,8,8,8,8};
	int PD_plan_size = 10;

	
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	
	shift = 0;
	output_shift = 0;
	startTaps = 0;
	iteration=0;
	
	max_iteration = Get_max_iteration(max_boxcar_width, PD_plan, PD_plan_size);
	
	// ----------> First iteration
	Get_next_iteration_parameters(0, PD_plan[0], nDMs, nTimesamples, &nBlocks, &unprocessed_samples, &decimated_timesamples, &shift, &output_shift, &startTaps, &iteration);
	nBoxcars=PD_plan[0];
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	PD_GPU_1st_float1<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars);
	
	// ----------> Higher iteration
	for(int f=1; f<PD_plan_size; f++){
		if(f<max_iteration){
			Get_next_iteration_parameters(nBoxcars, PD_plan[f], nDMs, nTimesamples, &nBlocks, &unprocessed_samples, &decimated_timesamples, &shift, &output_shift, &startTaps, &iteration);
			nBoxcars=PD_plan[f];
			gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
			blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
			if( (f%2) == 0 ) {
				PD_GPU_Nth_float1<<<gridSize,blockSize>>>(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[output_shift], &d_output_taps[output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration));
			}
			else {
				PD_GPU_Nth_float1<<<gridSize,blockSize>>>(&d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[output_shift], &d_output_taps[output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration));
			}
		}
	}

	*t_max_iterarion=max_iteration;
	return(unprocessed_samples);
}

