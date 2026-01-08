//Added by Karel Adamek
//#define HS_DEBUG

#include "aa_params.hpp"
#include "aa_device_harmonic_summing_kernel.hpp"
#include <stdio.h>

namespace astroaccelerate {

  /** \brief Computes a 1D simple harmonic sum that adds integer multiples of fundamental frequency only. */
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
    call_kernel_simple_harmonic_sum_GPU_kernel(gridSize, blockSize, d_input, d_output_SNR, d_output_harmonics, d_MSD, nTimesamples, nDMs, nHarmonics);
    
    return(0);
  }

  /** \brief Computes a 1D greedy harmonic sum that accounts for frequency shift in higher harmonics. */
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
    int nThreads = HRMS_ConstParams::nThreads;
    
    //---------> Task specific
    int nBlocks_x, nBlocks_y;
    nBlocks_x = (nTimesamples + nThreads - 1)/nThreads;
    nBlocks_y = nDMs;
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(nThreads, 1, 1);
    
    #ifdef HS_DEBUG
    if(DEBUG) printf("Data dimensions: %zu x %zu;\n",nDMs, nTimesamples);
    if(DEBUG) printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    if(DEBUG) printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
    #endif
    
    //---------> Greedy harmonic sum
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
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
  
  /** \brief Computes a 1D harmonic sum that is similar to the harmonic sum in PRESTO but adds more harmonics. */
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
    int nThreads = HRMS_ConstParams::nThreads;
    
    //---------> Task specific
    int nBlocks_x, nBlocks_y;
    nBlocks_x = (nTimesamples + nThreads - 1)/nThreads;
    nBlocks_y = nDMs;
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(nThreads, 1, 1);
    
    #ifdef HS_DEBUG
    if(DEBUG) printf("Data dimensions: %zu x %zu;\n",nDMs, nTimesamples);
    if(DEBUG) printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    if(DEBUG) printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
    #endif
    
    //---------> PRESTO harmonic sum
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
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
  
  /** \brief Computes a 1D harmonic sum that is similar to the harmonic sum in PRESTO. */
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
    int nThreads = HRMS_ConstParams::nThreads;
    
    //---------> Task specific
    int nBlocks_x, nBlocks_y;
    nBlocks_x = (nTimesamples + nThreads - 1)/nThreads;
    nBlocks_y = nDMs;
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(nThreads, 1, 1);
    int nHarmonicsFactor = (int) (log(nHarmonics)/log(2.0)) + 1;
    
    #ifdef HS_DEBUG
    if(DEBUG) printf("Data dimensions: %zu x %zu;\n",nDMs, nTimesamples);
    if(DEBUG) printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    if(DEBUG) printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
    #endif
    
    //---------> PRESTO harmonic sum
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    call_kernel_presto_harmonic_sum_GPU_kernel(
        gridSize,
        blockSize,
        d_input,
        d_output_SNR,
        d_output_harmonics,
        d_MSD,
        nTimesamples,
        nDMs,
        nHarmonicsFactor,
        enable_scalloping_loss_removal
    );
    
    return(0);
  }
  
  /** \brief Computes a 2D greedy harmonic sum. Not working properly at the moment. */
int periodicity_two_dimensional_greedy_harmonic_summing(
      float *d_input,
      float *d_output_max,
      float *d_output_SNR,
      ushort *d_output_harmonics,
      float *d_MSD,
      size_t N_f,
      size_t N_fdot,
      size_t max_f_idx,
      size_t max_fdot_idx,
      size_t nHarmonics
  ) {
    int nThreads = HRMS_ConstParams::nThreads;

    //---------> Task specific
    int nBlocks = (N_f * N_fdot + nThreads - 1) / nThreads;
    dim3 gridSize(nBlocks);
    dim3 blockSize(nThreads);

    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    call_two_dimensional_greedy_harmonic_sum_GPU_kernel(
        gridSize,
        blockSize,
        d_input,
        d_output_max,
        d_output_SNR,
        d_output_harmonics,
        d_MSD,
        N_f,
        N_fdot,
        max_f_idx,
        max_fdot_idx,
        nHarmonics
    );

    return (0);
  }
  

/**
 * @brief Host wrapper around @D greedy harmonic sum that configures the GPU kernel and executes it.
 *
 * Array dimensions are 2D f-fdot planes.
 *
 * @param d_output_pow Output 2D array of power values for highest SNR detected for each r and z.
 * @param d_output_SNR Output 2D array containing highest SNR value detected for each r and z.
 * @param d_output_harmonics Output 2D array containing number of harmonics summed to get returned SNR for each r and z.
 * @param d_output_shifts Output 2D array containing shifts in r and z to get returned SNR for each r and z.
 * @param d_input Input f-fdot plane
 * @param d_MSD Input mean and standard deviation.
 * @param num_r_bins Number of r bins.
 * @param num_z_bins Number of z bins.
 * @param nHarmonics Number of maximum number of harmonics to sum.
 */
int fdas_greedy_harmonic_summing_2d(
      float *d_output_power,
      float *d_output_SNR,
      short int *d_output_harmonics,
      short int *d_output_shifts,
      float *d_input,
      float *d_MSD,
      int num_r_bins,
      int num_z_bins,
      int nHarmonics
  ) {
    int nThreads = 128;
    int nBlocks_x = (num_r_bins + nThreads - 1) / nThreads;
    int nBlocks_y = num_z_bins;
    int nBlocks_z = 1;
    dim3 gridSize(nBlocks_x, nBlocks_y, nBlocks_z);
    dim3 blockSize(nThreads, 1, 1);

    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

    call_greedy_harmonic_sum_2d_kernel(
        gridSize,
        blockSize,
        d_output_power, 
        d_output_SNR, 
        d_output_harmonics, 
        d_output_shifts,
        d_input, 
        d_MSD, 
        num_r_bins, 
        num_z_bins, 
        nHarmonics
    );

    return (0);
  }
} //namespace astroaccelerate
