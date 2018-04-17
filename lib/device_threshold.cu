//Added by Karel Adamek
#define THRESHOLD_DEBUG

#include <helper_cuda.h>

#include "headers/device_BC_plan.h"

#include "headers/params.h"
#include "device_threshold_kernel.cu"

void THR_init(void) {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}

int THRESHOLD(float *d_input, ushort *d_input_taps, float *d_output_list, int *gmem_pos, float threshold, int nDMs, int nTimesamples, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int max_list_size, cudaStream_t streams) {
	//---------> Task specific
	int nBlocks, nRest, Elements_per_block, output_offset, decimated_timesamples, local_offset;
	int nCUDAblocks_x, nCUDAblocks_y;
	
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(WARP*THR_WARPS_PER_BLOCK, 1, 1);
	
	THR_init();
	
	output_offset=0;
	for(int f=0; f<max_iteration; f++){
		decimated_timesamples = PD_plan->operator[](f).decimated_timesamples;
		//local_offset = (offset>>f);
		local_offset = PD_plan->operator[](f).unprocessed_samples;
		if( (decimated_timesamples-local_offset)>0 ){
			Elements_per_block = WARP*THR_ELEM_PER_THREAD;
			nBlocks = (decimated_timesamples-local_offset)/Elements_per_block;
			nRest = (decimated_timesamples-local_offset) - nBlocks*Elements_per_block;
			if(nRest>0) nBlocks++;
			
			nCUDAblocks_x = nBlocks;
			nCUDAblocks_y = nDMs/THR_WARPS_PER_BLOCK;
			
			gridSize.x=nCUDAblocks_x; gridSize.y=nCUDAblocks_y; gridSize.z=1;
			blockSize.x=WARP*THR_WARPS_PER_BLOCK; blockSize.y=1; blockSize.z=1;
			
			output_offset = nDMs*PD_plan->operator[](f).output_shift;
			
			THR_GPU_WARP<<<gridSize, blockSize, 0, streams>>>(&d_input[output_offset], &d_input_taps[output_offset], d_output_list, gmem_pos, threshold, decimated_timesamples, decimated_timesamples-local_offset, shift, max_list_size, (1<<f));
			
			checkCudaErrors(cudaGetLastError());
		}
	}
	return (0);
}


int Threshold_for_periodicity_old(float *d_input, ushort *d_input_harms, float *d_output_list,  int *gmem_pos, float *d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int inBin, int max_list_size) {
	//---------> Task specific
	int nBlocks_p, nBlocks_s;
	
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(WARP, WARP/2, 1);
	
	nBlocks_p = (int) (primary_size/(blockSize.x*THR_ELEM_PER_THREAD));
	if( (primary_size%(blockSize.x*THR_ELEM_PER_THREAD))!=0 ) nBlocks_p++;
	
	nBlocks_s = (int) (secondary_size/blockSize.y);
	if( (secondary_size%blockSize.y)!=0 ) nBlocks_s++;
	
	gridSize.x = nBlocks_p;
	gridSize.y = nBlocks_s;
	gridSize.z = 1;
	
	#ifdef THRESHOLD_DEBUG
	printf("Primary:%d; Secondary:%d\n", primary_size, secondary_size);
	printf("gridSize: [%d; %d; %d]\n", gridSize.x, gridSize.y, gridSize.z);
	printf("blockSize: [%d; %d; %d]\n", blockSize.x, blockSize.y, blockSize.z);
	#endif
	
	THR_init();
	GPU_Threshold_for_periodicity_kernel_old<<<gridSize, blockSize>>>(d_input, d_input_harms, d_output_list, gmem_pos, d_MSD, threshold, primary_size, secondary_size, DM_shift, max_list_size, inBin);
	
	checkCudaErrors(cudaGetLastError());

	return (0);
}

int Threshold_for_periodicity(float *d_input, ushort *d_input_harms, float *d_output_list,  int *gmem_pos, float *d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int inBin, int max_list_size) {
	//---------> Task specific
	int nBlocks_p, nBlocks_s;
	
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(WARP, WARP/2, 1);
	
	nBlocks_p = (int) (primary_size/(blockSize.x*THR_ELEM_PER_THREAD));
	if( (primary_size%(blockSize.x*THR_ELEM_PER_THREAD))!=0 ) nBlocks_p++;
	
	nBlocks_s = (int) (secondary_size/blockSize.y);
	if( (secondary_size%blockSize.y)!=0 ) nBlocks_s++;
	
	gridSize.x = nBlocks_p;
	gridSize.y = nBlocks_s-1;
	gridSize.z = 1;
	
	#ifdef THRESHOLD_DEBUG
	printf("Primary:%d; Secondary:%d\n", primary_size, secondary_size);
	printf("gridSize: [%d; %d; %d]\n", gridSize.x, gridSize.y, gridSize.z);
	printf("blockSize: [%d; %d; %d]\n", blockSize.x, blockSize.y, blockSize.z);
	#endif
	
	THR_init();

	GPU_Threshold_for_periodicity_kernel<<<gridSize, blockSize>>>(d_input, d_input_harms, d_output_list, gmem_pos, d_MSD, threshold, primary_size, secondary_size, DM_shift, max_list_size, inBin);
	
	checkCudaErrors(cudaGetLastError());

	return (0);
}
