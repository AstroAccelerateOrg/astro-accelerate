//Added by Karel Adamek
//#define SPS_LONG_DEBUG

#include "aa_device_SPS_long.hpp"

namespace astroaccelerate {

  size_t Get_memory_requirement_of_SPS(){
    return((size_t) (5.5*sizeof(float) + 2*sizeof(ushort)));
  }


  void Assign_parameters(int f, std::vector<PulseDetection_plan> *PD_plan, int *decimated_timesamples, int *dtm, int *iteration, int *nBoxcars, int *nBlocks, int *output_shift, int *shift, int *startTaps, int *unprocessed_samples, int *total_ut){
    *decimated_timesamples = PD_plan->operator[](f).decimated_timesamples;
    *dtm                   = PD_plan->operator[](f).dtm;
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


  int SPDT_search_long_MSD_plane(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD_interpolated, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nTimesamples, int nDMs) {
    //---------> CUDA block and CUDA grid parameters
    dim3 gridSize(1, 1, 1);
    dim3 blockSize(PD_NTHREADS, 1, 1);
	
    //---------> Pulse detection FIR
    PD_SEARCH_LONG_init();
	
    int f;
    int decimated_timesamples, dtm, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut, MSD_plane_pos;
	
    // ----------> First iteration
    Assign_parameters(0, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
    MSD_plane_pos = 0;
    gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
    blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	
#ifdef SPS_LONG_DEBUG
    printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d; MSD_plane_pos:%d;\n",decimated_timesamples, dtm, iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut, MSD_plane_pos);
#endif
	
    if(nBlocks>0) call_kernel_SPDT_GPU_1st_plane(gridSize, blockSize, d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, (float2 *) d_MSD_interpolated, decimated_timesamples, nBoxcars, dtm);
	
    //checkCudaErrors(cudaGetLastError());
	
    for(f=1; f<max_iteration; f++){
      MSD_plane_pos = MSD_plane_pos + nBoxcars;
      Assign_parameters(f, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
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
		
      //checkCudaErrors(cudaGetLastError());
    }

    return(0);
  }

} //namespace astroaccelerate
