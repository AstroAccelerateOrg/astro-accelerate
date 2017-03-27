//Added by Karel Adamek

#include <vector>

#include "headers/params.h"
#include "headers/device_BC_plan.h"
#include "device_SPS_long_kernel.cu"


void Get_next_iteration_parameters_IF(int old_nBoxcars, int new_nBoxcars, int nDMs, int nTimesamples, int *nBlocks, int *unprocessed_samples, int *decimated_timesamples, int *shift, int *output_shift, int *startTaps, int  *iteration) {
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
	t_nBlocks=(t_decimated_timesamples)/t_Elements_per_block;
	t_nRest=t_decimated_timesamples - t_nBlocks*t_Elements_per_block;
	if(t_nRest>0) t_nBlocks++;
	
	*decimated_timesamples = t_decimated_timesamples;
	*unprocessed_samples = (*unprocessed_samples)/2 + new_nBoxcars + 6;
	*nBlocks = t_nBlocks;
}

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
	t_nBlocks=(t_decimated_timesamples - (*unprocessed_samples)/2 - new_nBoxcars)/t_Elements_per_block;
	t_nRest=t_decimated_timesamples - t_nBlocks*t_Elements_per_block;
	
	*decimated_timesamples = t_decimated_timesamples;
	*unprocessed_samples = t_nRest;
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

void Assign_parameters(int f, std::vector<PulseDetection_plan> *PD_plan, int *decimated_timesamples, int *iteration, int *nBoxcars, int *nBlocks, int *output_shift, int *shift, int *startTaps, int *unprocessed_samples, int *total_ut){
	*decimated_timesamples = PD_plan->operator[](f).decimated_timesamples;
	*iteration             = PD_plan->operator[](f).iteration;
	*nBoxcars              = PD_plan->operator[](f).nBoxcars;
	*nBlocks               = PD_plan->operator[](f).nBlocks;
	*output_shift          = PD_plan->operator[](f).output_shift;
	*shift                 = PD_plan->operator[](f).shift;           
	*startTaps             = PD_plan->operator[](f).startTaps; 
	*unprocessed_samples   = PD_plan->operator[](f).unprocessed_samples;
	*total_ut              = PD_plan->operator[](f).total_ut;	
}

void PD_SEARCH_LONG_init() {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeEightByte);
}

int PD_SEARCH_LONG(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int max_boxcar_width, int nDMs, int nTimesamples, int *t_max_iterarion) {
	//---------> Task specific
	int nBlocks, unprocessed_samples;//nRest, Elements_per_block,
	int shift, output_shift, iteration, max_iteration, startTaps, decimated_timesamples;
	int nBoxcars;
	int PD_plan[10]={32,16,16,16,8,8,8,8,8,8};
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

int PD_SEARCH_LONG_BLN(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int max_boxcar_width, int nDMs, int nTimesamples, int *t_max_iterarion) {
	//---------> Task specific
	int nBlocks, unprocessed_samples;//nRest, Elements_per_block,
	int shift, output_shift, iteration, max_iteration, startTaps, decimated_timesamples;
	int nBoxcars;
	int PD_plan[10]={32,16,16,16,8,8,8,8,8,8};
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
	printf("decimated_timesamples:%d; iteration:%d; nBlocks:%d; nBoxcars:%d; shift:%d; start_Taps:%d;\n", decimated_timesamples, iteration, gridSize.x, nBoxcars, shift, startTaps);
	PD_GPU_1st_float1_BLN<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars);
	
	// ----------> Higher iteration
	for(int f=1; f<PD_plan_size; f++){
		if(f<max_iteration){
			Get_next_iteration_parameters(nBoxcars, PD_plan[f], nDMs, nTimesamples, &nBlocks, &unprocessed_samples, &decimated_timesamples, &shift, &output_shift, &startTaps, &iteration);
			nBoxcars=PD_plan[f];
			gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
			blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
			printf("decimated_timesamples:%d; iteration:%d; nBlocks:%d; nBoxcars:%d; shift:%d; start_Taps:%d;\n", decimated_timesamples, iteration, gridSize.x, nBoxcars, shift, startTaps);
			if( (f%2) == 0 ) {
				PD_GPU_Nth_float1_BLN<<<gridSize,blockSize>>>(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[output_shift], &d_output_taps[output_shift], decimated_timesamples, nBoxcars, startTaps, (1<<iteration));
			}
			else {
				PD_GPU_Nth_float1_BLN<<<gridSize,blockSize>>>(&d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[output_shift], &d_output_taps[output_shift], decimated_timesamples, nBoxcars, startTaps, (1<<iteration));
			}
		}
	}

	*t_max_iterarion=max_iteration;
	unprocessed_samples = unprocessed_samples*(1<<(max_iteration-1));
	return(unprocessed_samples);
}


 int PD_SEARCH_LONG_BLN_IF(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples) {
	//---------> Task specific
	
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	int f;
	int decimated_timesamples, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut;
	
	// ----------> First iteration
	Assign_parameters(0, PD_plan, &decimated_timesamples, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	printf("decimated_timesamples:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples	,iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
	PD_GPU_1st_float1_BLN_IF<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars);
	
	
	for(f=1; f<max_iteration; f++){
		Assign_parameters(f, PD_plan, &decimated_timesamples, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
		gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
		blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
		printf("decimated_timesamples:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples	,iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
		if( (f%2) == 0 ) {
			PD_GPU_Nth_float1_BLN_IF<<<gridSize,blockSize>>>(d_input, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), shift);
		}
		else {
			PD_GPU_Nth_float1_BLN_IF<<<gridSize,blockSize>>>(d_decimated, d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), shift);
		}
	}
	
	
	/*
	int nBlocks, unprocessed_samples;
	int shift, output_shift, iteration, startTaps, decimated_timesamples;
	int nBoxcars;
	
	shift = 0;
	output_shift = 0;
	startTaps = 0;
	iteration=0;
	
	// ----------> First iteration
	Get_next_iteration_parameters_IF(0, BC_widths[0], nDMs, nTimesamples, &nBlocks, &unprocessed_samples, &decimated_timesamples, &shift, &output_shift, &startTaps, &iteration);
	nBoxcars=BC_widths[0];
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	printf("decimated_timesamples:%d; iteration:%d; nBlocks:%d; nBoxcars:%d; shift:%d; start_Taps:%d;\n", decimated_timesamples, iteration, gridSize.x, nBoxcars, shift, startTaps);
	PD_GPU_1st_float1_BLN_IF<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars);
	
	
	for(int f=1; f< (int) BC_widths_size; f++){
		if(f<max_iteration){
			Get_next_iteration_parameters_IF(nBoxcars, BC_widths[f], nDMs, nTimesamples, &nBlocks, &unprocessed_samples, &decimated_timesamples, &shift, &output_shift, &startTaps, &iteration);
			nBoxcars=BC_widths[f];
			gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
			blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
			printf("decimated_timesamples:%d; iteration:%d; nBlocks:%d; nBoxcars:%d; shift:%d; start_Taps:%d;\n", decimated_timesamples, iteration, gridSize.x, nBoxcars, shift, startTaps);
			if( (f%2) == 0 ) {
				PD_GPU_Nth_float1_BLN_IF<<<gridSize,blockSize>>>(d_input, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[output_shift], &d_output_taps[output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), shift);
			}
			else {
				PD_GPU_Nth_float1_BLN_IF<<<gridSize,blockSize>>>(d_decimated, d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[output_shift], &d_output_taps[output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), shift);
			}
		}
	}
	*/

	unprocessed_samples = unprocessed_samples*(1<<(max_iteration-1));
	return(unprocessed_samples);
}




