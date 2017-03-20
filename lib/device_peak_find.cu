// James Sharpe's peak finding code

#include "headers/device_BC_plan.h"

#include <npp.h>
#include "helper_cuda.h"

#include "headers/params.h"
#include "device_peak_find_kernel.cu"


void PEAK_FIND(float *d_output_SNR, ushort *d_output_taps, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration){
	int output_offset, decimated_timesamples, local_offset;

        cudaStream_t stream;
	dim3 blockDim(32, 2, 1);
	dim3 gridSize(1 + ((nTimesamples-1)/blockDim.x), 1 + ((nDMs-1)/blockDim.y), 1);
 
        //Queue up each iteration on a stream so that they execute with minimal time between them
        //Note in theory each call could be on a separate stream to allow concurrent kernel execution
        //however this is making an assumption that atomicAdd on the same memory from two concurrent kernels is ok
        // and this kernel is limited by memory bandwidth so running concurrently is not that important
        cudaStreamCreate(&stream);
	
	output_offset=0;
	for(int f=0; f<max_iteration; f++){
		decimated_timesamples = (nTimesamples>>f);
		//local_offset = offset>>f;
		local_offset = PD_plan->operator[](f).unprocessed_samples;
		if( (decimated_timesamples-local_offset-1) > 0 ) {
			gridSize.x = 1 + ((decimated_timesamples-local_offset-1)/blockDim.x);
			gridSize.y = 1 + ((nDMs-1)/blockDim.y);
			gridSize.z = 1;
			
			dilate_peak_find<<<gridSize, blockDim, 0, stream>>>(&d_output_SNR[output_offset], &d_output_taps[output_offset], (float4*)d_peak_list, decimated_timesamples, nDMs, local_offset, threshold, max_peak_size, gmem_peak_pos, shift, (1<<f));
			output_offset = output_offset + nDMs*decimated_timesamples;
		}
	}

        cudaStreamSynchronize(stream);
        cudaStreamDestroy(stream);
}

void PEAK_FIND_FOR_FDAS(float *d_ffdot_plane, float *d_peak_list, float *d_MSD, int nDMs, int nTimesamples, float threshold, unsigned int max_peak_size, unsigned int *gmem_peak_pos, float DM_trial) {
	dim3 blockDim(32, 2, 1);
	dim3 gridSize(1 + ((nTimesamples-1)/blockDim.x), 1 + ((nDMs-1)/blockDim.y), 1);
	
	dilate_peak_find_for_fdas<<<gridSize, blockDim>>>( d_ffdot_plane, d_peak_list, d_MSD, nTimesamples, nDMs, 0, threshold, max_peak_size, gmem_peak_pos, DM_trial);
}