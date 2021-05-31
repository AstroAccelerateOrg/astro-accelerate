//Added by Karel Adamek
//#define HS_DEBUG

#include "aa_params.hpp"
#include "aa_device_harmonic_summing_kernel.hpp"

namespace astroaccelerate {

  /** \brief Computes a simple harmonic sum for periodicity. */
  int periodicity_simple_harmonic_summing(
      float *d_input, 
      float *d_output_SNR,
      ushort *d_output_harmonics,
      float *d_MSD,
      int nTimesamples,
      int nDMs,
      int nHarmonics
  ){
    //---------> Task specific
    int nBlocks_x, nBlocks_y;
    nBlocks_x = nTimesamples;
    nBlocks_y = (nDMs + PHS_NTHREADS - 1)/PHS_NTHREADS;
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(PHS_NTHREADS, 1, 1);
    
    #ifdef HS_DEBUG
    printf("Data dimensions: %d x %d;\n", nDMs, nTimesamples);
    printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
    #endif
    
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    call_kernel_simple_harmonic_sum_GPU_kernel(gridSize, blockSize, d_input, d_output_SNR, d_output_harmonics, d_MSD, nTimesamples, nDMs, nHarmonics);
    
    return(0);
  }

  int periodicity_greedy_harmonic_summing(
      float *d_input,
      float *d_output_SNR,
      ushort *d_output_harmonics,
      float *d_MSD,
      int nTimesamples,
      int nDMs,
      int nHarmonics,
      int enable_scalloping_loss_removal
  ){
    //---------> Task specific
    int nBlocks_x, nBlocks_y;
    nBlocks_x = (nTimesamples + GHRMS_NTHREADS - 1)/GHRMS_NTHREADS;
    nBlocks_y = nDMs;
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(GHRMS_NTHREADS, 1, 1);
    
    #ifdef HS_DEBUG
    if(DEBUG) printf("Data dimensions: %zu x %zu;\n",nDMs, nTimesamples);
    if(DEBUG) printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    if(DEBUG) printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
    #endif
    
    //---------> Greedy harmonic sum
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    call_kernel_greedy_harmonic_sum_GPU_kernel(
        gridSize,
        blockSize,
        d_input,
        d_output_SNR,
        d_output_harmonics,
        d_MSD,
        nTimesamples,
        nDMs,
        nHarmonics,
        enable_scalloping_loss_removal
    );
    
    return(0);
  }
  
  int periodicity_presto_plus_harmonic_summing(
      float *d_input,
      float *d_output_SNR,
      ushort *d_output_harmonics,
      float *d_MSD,
      int nTimesamples,
      int nDMs,
      int nHarmonics,
      int enable_scalloping_loss_removal
  ) {
    //---------> Task specific
    int nBlocks_x, nBlocks_y;
    nBlocks_x = (nTimesamples + GHRMS_NTHREADS - 1)/GHRMS_NTHREADS;
    nBlocks_y = nDMs;
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(GHRMS_NTHREADS, 1, 1);
    
    #ifdef HS_DEBUG
    if(DEBUG) printf("Data dimensions: %zu x %zu;\n",nDMs, nTimesamples);
    if(DEBUG) printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    if(DEBUG) printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
    #endif
    
    //---------> PRESTO harmonic sum
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    call_kernel_presto_plus_harmonic_sum_GPU_kernel(
        gridSize,
        blockSize,
        d_input,
        d_output_SNR,
        d_output_harmonics,
        d_MSD,
        nTimesamples,
        nDMs,
        nHarmonics,
        enable_scalloping_loss_removal
    );
    
    return(0);
  }
  
  /*
  int periodicity_presto_harmonic_summing(
      float *d_input,
      float *d_output_SNR,
      ushort *d_output_harmonics,
      float *d_MSD,
      int nTimesamples,
      int nDMs,
      int nHarmonics,
      int enable_scalloping_loss_removal
  ) {
    //---------> Task specific
    int nBlocks_x, nBlocks_y;
    nBlocks_x = (nTimesamples + GHRMS_NTHREADS - 1)/GHRMS_NTHREADS;
    nBlocks_y = nDMs;
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(GHRMS_NTHREADS, 1, 1);
    
    #ifdef HS_DEBUG
    if(DEBUG) printf("Data dimensions: %zu x %zu;\n",nDMs, nTimesamples);
    if(DEBUG) printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    if(DEBUG) printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
    #endif
    
    //---------> PRESTO harmonic sum
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    call_kernel_presto_harmonic_sum_GPU_kernel(
        gridSize,
        blockSize,
        d_input,
        d_output_SNR,
        d_output_harmonics,
        d_MSD,
        nTimesamples,
        nDMs,
        nHarmonics,
        enable_scalloping_loss_removal
    );
    
    return(0);
  }
  */
  
} //namespace astroaccelerate
