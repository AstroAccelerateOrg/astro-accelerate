//Added by Karel Adamek
#include "headers/device_BC_plan.h"

#include "headers/params.h"
#include "device_threshold_kernel.cu"

void THR_init(void) {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}

int THRESHOLD(float *d_input, ushort *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nDMs, int nTimesamples, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int max_list_size) {
	//---------> Task specific
	int nBlocks, nRest, Elements_per_block, output_offset, decimated_timesamples, local_offset;
	int nCUDAblocks_x, nCUDAblocks_y;
	
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(WARP*THR_WARPS_PER_BLOCK, 1, 1);
	
	THR_init();
	
	output_offset=0;
	for(int f=0; f<max_iteration; f++){
		decimated_timesamples = (nTimesamples>>f);
		//local_offset = (offset>>f);
		local_offset = PD_plan->operator[](f).unprocessed_samples;
		if( (decimated_timesamples-local_offset)>0 ){
			Elements_per_block = 2*WARP*THR_ELEM_PER_THREAD;
			nBlocks = (decimated_timesamples-local_offset)/Elements_per_block;
			nRest = (decimated_timesamples-local_offset) - nBlocks*Elements_per_block;
			if(nRest>0) nBlocks++;
			
			nCUDAblocks_x = nBlocks;
			nCUDAblocks_y = nDMs/THR_WARPS_PER_BLOCK;
			
			gridSize.x=nCUDAblocks_x; gridSize.y=nCUDAblocks_y; gridSize.z=1;
			blockSize.x=WARP*THR_WARPS_PER_BLOCK; blockSize.y=1; blockSize.z=1;
			
			THR_GPU_WARP<<<gridSize, blockSize>>>((float2 *) &d_input[output_offset], &d_input_taps[output_offset], d_output_list, gmem_pos, threshold, decimated_timesamples, decimated_timesamples-local_offset, shift, max_list_size, (1<<f));
			output_offset = output_offset + nDMs*decimated_timesamples;
		}
	}

	return (0);
}
