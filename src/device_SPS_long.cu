//Added by Karel Adamek
//#define SPS_LONG_DEBUG

#include "device_SPS_long.hpp"
#include "device_SPS_Plan.hpp"


size_t Get_memory_requirement_of_SPS(){
	return((size_t) (5.5*sizeof(float) + 2*sizeof(ushort)));
}


void Assign_parameters(ProcessingDetails &details, int *decimated_timesamples, int *dtm, int *iteration, int *nBoxcars, int *nBlocks, int *output_shift, int *shift, int *startTaps, int *unprocessed_samples, int *total_ut) {
	*decimated_timesamples = details.decimated_timesamples;
	*dtm                   = details.dtm;
	*iteration             = details.iteration;
	*nBlocks               = details.number_blocks;
	*nBoxcars              = details.number_boxcars;
	*output_shift          = details.output_shift;
	*shift                 = details.shift;           
	*startTaps             = details.start_taps; 
	*total_ut              = details.total_unprocessed;
	*unprocessed_samples   = details.unprocessed_samples;
}

void PD_SEARCH_LONG_init() {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeEightByte);
}


int SPDT_search_long_MSD_plane(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD_interpolated, SPS_Plan &spsplan, int max_iteration, int nTimesamples, int nDMs) {
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	int f;
	int decimated_timesamples, dtm, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut, MSD_plane_pos;
	
	// ----------> First iteration
	Assign_parameters(spsplan.GetDetails(0), &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
	MSD_plane_pos = 0;
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	
	#ifdef SPS_LONG_DEBUG
	printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d; MSD_plane_pos:%d;\n",decimated_timesamples, dtm, iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut, MSD_plane_pos);
	#endif
	
	if(nBlocks>0) call_kernel_SPDT_GPU_1st_plane(gridSize, blockSize, d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, (float2 *) d_MSD_interpolated, decimated_timesamples, nBoxcars, dtm);
	
	checkCudaErrors(cudaGetLastError());
	
	for(f=1; f<max_iteration; f++){
		MSD_plane_pos = MSD_plane_pos + nBoxcars;
		// TODO: Need to update the spsplan here
		Assign_parameters(spsplan.GetDetails(f), &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
		gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
		blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
		
		#ifdef SPS_LONG_DEBUG
		printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d; MSD_plane_pos:%d;\n",decimated_timesamples, dtm, iteration, nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut, MSD_plane_pos);
		#endif
		
		if( (f%2) == 0 ) {
			if(nBlocks>0) 
			  call_kernel_SPDT_GPU_Nth_plane(gridSize,blockSize, &d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], (float2 *) &d_MSD_interpolated[MSD_plane_pos*2], decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		else {
			if(nBlocks>0) 
			  call_kernel_SPDT_GPU_Nth_plane(gridSize,blockSize, &d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], (float2 *) &d_MSD_interpolated[MSD_plane_pos*2], decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		
		checkCudaErrors(cudaGetLastError());
	}

	return(0);
}

