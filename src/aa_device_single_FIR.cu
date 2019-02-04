//Added by Karel Adamek

#include "aa_params.hpp"
#include "aa_device_single_FIR_kernel.hpp"

namespace astroaccelerate {

  void PD_FIR_init(void)
  {
    //---------> Specific nVidia stuff
    cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
    cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeEightByte);
  }

  int PD_FIR(float *d_input, float *d_output, int nTaps, int nDMs, int nTimesamples)
  {
    //---------> Task specific
    int ut; //unused timesamples
    int itemp = (int) ( ( nTaps - 1 ) / ( WARP * PD_FIR_ACTIVE_WARPS ) ) + 1;
    int nLoops = PD_FIR_NWINDOWS + itemp;

    //---------> CUDA block and CUDA grid parameters
    int nCUDAblocks_x = (int) ( ( nTimesamples - nTaps + 1 ) / ( PD_FIR_ACTIVE_WARPS * WARP * PD_FIR_NWINDOWS ) );
    int nCUDAblocks_y = nDMs;
    int SM_size = ( PD_FIR_ACTIVE_WARPS * WARP * PD_FIR_NWINDOWS + nTaps - 1 ) * 4;

    dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);
    dim3 blockSize(PD_FIR_ACTIVE_WARPS * WARP, 1, 1);

    //---------> Pulse detection FIR
    PD_FIR_init();
    call_kernel_PD_FIR_GPU(gridSize, blockSize, SM_size, d_input, d_output, nTaps, nLoops, nTimesamples);

    ut = nTimesamples - nCUDAblocks_x * PD_FIR_ACTIVE_WARPS * WARP * PD_FIR_NWINDOWS;
    return ( ut );
  }


  int GPU_FIRv1_wrapper(float *d_input, float *d_output, int nTaps, unsigned int nDMs, unsigned int nTimesamples){
    //---------> Task specific
    int ut; //unused timesamples
    int itemp=(int) ((nTaps - 1)/(WARP*PD_FIR_ACTIVE_WARPS)) + 1;
    int nLoops=PD_FIR_NWINDOWS + itemp;
	
    //---------> CUDA block and CUDA grid parameters
    int nCUDAblocks_x=(int) ((nTimesamples - nTaps + 1)/(PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS));
    int nCUDAblocks_y=nDMs; //Head size
    int SM_size=(PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS + nTaps - 1)*4;
	
    dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);			//nCUDAblocks_y goes through spectra
    dim3 blockSize(PD_FIR_ACTIVE_WARPS*WARP, 1, 1); 		//nCUDAblocks_x goes through channels
	
    // ----------------------------------------------->
    // --------> Measured part (Pulse detection FIR)	
	
    //---------> Pulse detection FIR
    cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    call_kernel_PD_FIR_GPUv1(gridSize,blockSize, SM_size, d_input, d_output, nTaps, nLoops, nTimesamples);

    // --------> Measured part (Pulse detection FIR)
    // ----------------------------------------------->
	
    ut=nTimesamples - nCUDAblocks_x*PD_FIR_ACTIVE_WARPS*WARP*PD_FIR_NWINDOWS;
    return(ut);
  }



  int PPF_L1(float *d_input, float *d_output, int nDMs, int nTimesamples, int nTaps) {
    //---------> CUDA block and CUDA grid parameters
    int itemp = nTimesamples-nTaps+1;
    int nCUDAblocks_x=(int) (itemp/PPF_L1_THREADS_PER_BLOCK);
    if(itemp%PPF_L1_THREADS_PER_BLOCK!=0) nCUDAblocks_x++;
    int nCUDAblocks_y=(int) nDMs;
    dim3 GridSize(nCUDAblocks_x, nCUDAblocks_y, 1);
    dim3 BlockSize(PPF_L1_THREADS_PER_BLOCK, 1, 1);
	
    //printf("nTimesamples:%d; nDMs:%d; nTaps:%d\n", nTimesamples, nDMs, nTaps);
    //printf("GridSize: [%d;%d;%d]\n", GridSize.x, GridSize.y, GridSize.z);
    //printf("BlockSize: [%d;%d;%d]\n", BlockSize.x, BlockSize.y, BlockSize.z);
	
    if( itemp>0 && GridSize.x>0 ){
      //printf("Running...\n");
      call_kernel_Fir_L1(GridSize, BlockSize, d_input, d_output, nTaps, nTimesamples);
    }
    else return(-1);
	
    return (nTaps-1);
  }

} //namespace astroaccelerate
