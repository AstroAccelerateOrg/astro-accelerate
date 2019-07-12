//Added by Karel Adamek
//#define THRESHOLD_DEBUG

#include <stdio.h>
#include "aa_device_threshold.hpp"
#include "aa_device_BC_plan.hpp"

#include "aa_params.hpp"
#include "aa_device_threshold_kernel.hpp"

namespace astroaccelerate {

  void THR_init(void) {
    //---------> Specific nVidia stuff
    cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
    cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
  }

  int SPDT_threshold(float *d_input, ushort *d_input_taps, unsigned int *d_output_list_DM, unsigned int *d_output_list_TS, float *d_output_list_SNR, unsigned int *d_output_list_BW, int *gmem_pos, float threshold, int nDMs, int nTimesamples, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int max_list_size) {
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
			
	call_kernel_THR_GPU_WARP(gridSize, blockSize, &d_input[output_offset], &d_input_taps[output_offset], d_output_list_DM, d_output_list_TS, d_output_list_SNR, d_output_list_BW, gmem_pos, threshold, decimated_timesamples, decimated_timesamples-local_offset, shift, max_list_size, (1<<f));
			
	//checkCudaErrors(cudaGetLastError());
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
    call_kernel_GPU_Threshold_for_periodicity_kernel_old(gridSize, blockSize, d_input, d_input_harms, d_output_list, gmem_pos, d_MSD, threshold, primary_size, secondary_size, DM_shift, max_list_size, inBin);
	
    //checkCudaErrors(cudaGetLastError());

    return (0);
  }

  int Threshold_for_periodicity(float *d_input, ushort *d_input_harms, float *d_output_list,  int *gmem_pos, float *d_MSD, float threshold, int primary_size, int secondary_size, int DM_shift, int inBin, int max_list_size) {
    //---------> Nvidia stuff
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, CARD);
    size_t max_x = deviceProp.maxGridSize[0], max_y = deviceProp.maxGridSize[1];
    THR_init();
	
    //---------> Task specific
    int nBlocks_p, nBlocks_s, nRepeats;
    size_t sec_per_chunk, sec_tail;
	
    dim3 gridSize(1, 1, 1);
    dim3 blockSize(WARP, WARP, 1);
	
    nBlocks_p = (int) (primary_size/(blockSize.x*THR_ELEM_PER_THREAD));
    if( (primary_size%(blockSize.x*THR_ELEM_PER_THREAD))!=0 ) nBlocks_p++;
    if((size_t) nBlocks_p > max_x) {printf("Too many DM trials!\n"); exit(1);}
	
    nBlocks_s = (int) (secondary_size/blockSize.y);
    if( (secondary_size%blockSize.y)!=0 ) nBlocks_s++;
	
    nRepeats = ceil( (float) ((float) nBlocks_s)/((float) max_y));
	
    sec_per_chunk = (int) (secondary_size/nRepeats);
    sec_tail = secondary_size - sec_per_chunk*(nRepeats-1);
	
    std::vector<size_t> secondary_size_per_chunk;
    for(int f=0; f<(nRepeats-1); f++){
      secondary_size_per_chunk.push_back(sec_per_chunk);
    }
    secondary_size_per_chunk.push_back(sec_tail);
	
#ifdef THRESHOLD_DEBUG
    printf("Primary:%d; Secondary:%d\n", primary_size, secondary_size);
    printf("gridSize: [%d; %d; %d]\n", nBlocks_p, nBlocks_s, gridSize.z);
    printf("blockSize: [%d; %d; %d]\n", blockSize.x, blockSize.y, blockSize.z);
    printf("nRepeats: %d; sec_per_chunk: %zu; sec_tail: %zu;\n", nRepeats, sec_per_chunk, sec_tail);
    printf("Secondary dimensions per chunk: ");
    for(size_t f=0; f<secondary_size_per_chunk.size(); f++) printf("%zu ", secondary_size_per_chunk[f]); printf("\n");	
#endif

    size_t shift = 0;
    for(int f=0; f<(int) secondary_size_per_chunk.size(); f++){
      nBlocks_s = (int) (secondary_size_per_chunk[f]/blockSize.y);
      if( (secondary_size_per_chunk[f]%blockSize.y)!=0 ) nBlocks_s++;
	
      gridSize.x = nBlocks_p;
      gridSize.y = nBlocks_s;
      gridSize.z = 1;
	
      call_kernel_GPU_Threshold_for_periodicity_kernel(gridSize, blockSize, &d_input[shift*primary_size], d_input_harms, d_output_list, gmem_pos, d_MSD, threshold, primary_size, secondary_size_per_chunk[f], DM_shift, max_list_size, inBin);
	
      //checkCudaErrors(cudaGetLastError());
      shift = shift + secondary_size_per_chunk[f];
    }

    return (0);
  }

} //namespace astroaccelerate
