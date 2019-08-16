// James Sharpe's peak finding code
//#define PEAK_FIND_DEBUG

#include <stdio.h>
#include "aa_device_peak_find.hpp"

namespace astroaccelerate {
  
	void SPDT_peak_find_stencil_7x7(float *d_output_SNR, ushort *d_output_taps, unsigned int *d_peak_list_DM, unsigned int *d_peak_list_TS, float *d_peak_list_SNR, unsigned int *d_peak_list_BW, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration){

		int decimated_timesamples, local_offset;

		for(int f=0; f<max_iteration; f++){
			decimated_timesamples = PD_plan->operator[](f).decimated_timesamples;
			local_offset = PD_plan->operator[](f).unprocessed_samples;

			if( (decimated_timesamples-local_offset-1)>0 ){
				int nThreads = 128;
				int nBlocks = decimated_timesamples*nDMs/nThreads + 1;

				call_kernel_peak_find_list(nBlocks, nThreads, &d_output_SNR[nDMs*PD_plan->operator[](f).output_shift], decimated_timesamples, nDMs, threshold, gmem_peak_pos, shift, (1<<f), &d_output_taps[nDMs*PD_plan->operator[](f).output_shift], max_peak_size, d_peak_list_DM, d_peak_list_TS, d_peak_list_SNR, d_peak_list_BW);
			} // if (decimated_timesamples ...)
		} //for max_iteration
	}

  void SPDT_peak_find(float *d_output_SNR, ushort *d_output_taps, unsigned int *d_peak_list_DM, unsigned int *d_peak_list_TS, float *d_peak_list_SNR, unsigned int *d_peak_list_BW, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration){
    int decimated_timesamples, local_offset;
	
    dim3 blockDim(32, 2, 1);
    dim3 gridSize(1, 1, 1);
	
    for(int f=0; f<max_iteration; f++){
      decimated_timesamples = PD_plan->operator[](f).decimated_timesamples;
      //local_offset = offset>>f;
      local_offset = PD_plan->operator[](f).unprocessed_samples;
      if( (decimated_timesamples-local_offset-1)>0 ){		
	gridSize.x = 1 + ((decimated_timesamples-local_offset-1)/blockDim.x);
	gridSize.y = 1 + ((nDMs-1)/blockDim.y);
	gridSize.z = 1;		
	
	call_kernel_dilate_peak_find(gridSize, blockDim, &d_output_SNR[nDMs*PD_plan->operator[](f).output_shift], &d_output_taps[nDMs*PD_plan->operator[](f).output_shift], d_peak_list_DM, d_peak_list_TS, d_peak_list_SNR, d_peak_list_BW, decimated_timesamples, nDMs, local_offset, threshold, max_peak_size, gmem_peak_pos, shift, (1<<f));
      }
    }
  }

  void PEAK_FIND_FOR_FDAS(float *d_ffdot_plane, float *d_peak_list, float *d_MSD, int nDMs, int nTimesamples, float threshold, unsigned int max_peak_size, unsigned int *gmem_peak_pos, float DM_trial) {
    dim3 blockDim(32, 2, 1);
    dim3 gridSize(1 + ((nTimesamples-1)/blockDim.x), 1 + ((nDMs-1)/blockDim.y), 1);
	
    call_kernel_dilate_peak_find_for_fdas(gridSize, blockDim, d_ffdot_plane, d_peak_list, d_MSD, nTimesamples, nDMs, 0, threshold, max_peak_size, gmem_peak_pos, DM_trial);
  }

  void Peak_find_for_periodicity_search_old(float *d_input_SNR, ushort *d_input_harmonics, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, float *d_MSD, int DM_shift, int inBin){
    dim3 blockDim(32, 32, 1);
    dim3 gridSize(1 + ((nTimesamples-1)/blockDim.x), 1 + ((nDMs-1)/blockDim.y), 1);
	
    call_kernel_dilate_peak_find_for_periods_old(gridSize, blockDim, d_input_SNR, d_input_harmonics, d_peak_list, nTimesamples, nDMs, 0, threshold, max_peak_size, gmem_peak_pos, d_MSD, DM_shift, inBin);
  }

  void Peak_find_for_periodicity_search(float *d_input_SNR, ushort *d_input_harmonics, float *d_peak_list, int secondary_size, int primary_size, float threshold, int max_peak_size, int *gmem_peak_pos, float *d_MSD, int DM_shift, int inBin){
    // nDMs = secondary_size
    // nTimesamples = primary_size
    //---------> Nvidia stuff
    // find maximum values for maximum grid sizes in x and y
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, CARD);
    size_t max_x = deviceProp.maxGridSize[0], max_y = deviceProp.maxGridSize[1];
	
    //---------> Task specific
    int nBlocks_p, nBlocks_s, nRepeats;
    size_t sec_per_chunk, sec_tail;
    dim3 blockDim(32, 32, 1);
	
    // number of thread-blocks required for the task
    nBlocks_p = 1 + ((primary_size-1)/blockDim.x);
    if((size_t) nBlocks_p > max_x) {printf("Too many DM trials!\n"); exit(1);}
    nBlocks_s = 1 + ((secondary_size-1)/blockDim.y);
	
    // number of kernel launches required for the task
    nRepeats = ceil( (float) ((float) nBlocks_s)/((float) max_y));
    sec_per_chunk = (int) (secondary_size/nRepeats);
    sec_tail = secondary_size - sec_per_chunk*(nRepeats-1);
	
    // creating vector with y dim processed by each kernel launch
    std::vector<size_t> secondary_size_per_chunk;
    for(int f=0; f<(nRepeats-1); f++){
      secondary_size_per_chunk.push_back(sec_per_chunk);
    }
    secondary_size_per_chunk.push_back(sec_tail);
	
#ifdef PEAK_FIND_DEBUG
    printf("Primary:%d; Secondary:%d\n", primary_size, secondary_size);
    printf("gridSize: [%d; %d; %d]\n", nBlocks_p, nBlocks_s, 1);
    printf("blockSize: [%d; %d; %d]\n", blockDim.x, blockDim.y, blockDim.z);
    printf("nRepeats: %d; sec_per_chunk: %zu; sec_tail: %zu;\n", nRepeats, sec_per_chunk, sec_tail);
    printf("Secondary dimensions per chunk: ");
    for(size_t f=0; f<secondary_size_per_chunk.size(); f++) printf("%zu ", secondary_size_per_chunk[f]); printf("\n");	
#endif
	
    // launching GPU kernels
    size_t shift = 0;
    for(int f=0; f<(int) secondary_size_per_chunk.size(); f++){
      nBlocks_s = 1 + ((secondary_size_per_chunk[f]-1)/blockDim.y);
      dim3 gridSize(nBlocks_p, nBlocks_s, 1);
	
      call_kernel_dilate_peak_find_for_periods(gridSize, blockDim, &d_input_SNR[shift*primary_size], d_input_harmonics, d_peak_list, primary_size, secondary_size, 0, threshold, max_peak_size, gmem_peak_pos, d_MSD, DM_shift, inBin);
		
      shift = shift + secondary_size_per_chunk[f];
    }
	
  }

} //namespace astroaccelerate
