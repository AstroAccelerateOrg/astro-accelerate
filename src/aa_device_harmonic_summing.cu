//Added by Karel Adamek
//#define HS_DEBUG

#include "aa_params.hpp"
#include "aa_device_harmonic_summing_kernel.hpp"

namespace astroaccelerate {

  /** \brief Computes a simple harmonic sum for periodicity. */
  void periodicity_simple_harmonic_summing_old(float *d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics){
    //---------> Task specific
    int nBlocks_x, nBlocks_y;

    nBlocks_x = nTimesamples;
    nBlocks_y = nSpectra/PHS_NTHREADS;
    if ( nSpectra%PHS_NTHREADS !=0 ) nBlocks_y++;
	
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(PHS_NTHREADS, 1, 1);
	
#ifdef HS_DEBUG
    printf("Data dimensions: %d x %d;\n", nSpectra, nTimesamples);
    printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
#endif
	
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    call_kernel_PHS_GPU_kernel_old(gridSize,blockSize, d_input, d_output_SNR, d_output_harmonics, d_MSD, nTimesamples, nSpectra, nHarmonics);
  }

  /** \brief Computes a simple harmonic sum for periodicity. */
  void periodicity_simple_harmonic_summing(float *d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics){
    //---------> Task specific
    int nBlocks_x, nBlocks_y;

    nBlocks_x = nTimesamples;
    nBlocks_y = nSpectra/PHS_NTHREADS;
    if ( nSpectra%PHS_NTHREADS !=0 ) nBlocks_y++;
	
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(PHS_NTHREADS, 1, 1);
	
#ifdef HS_DEBUG
    printf("Data dimensions: %d x %d;\n", nSpectra, nTimesamples);
    printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
#endif
	
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    call_kernel_PHS_GPU_kernel(gridSize, blockSize, d_input, d_output_SNR, d_output_harmonics, d_MSD, nTimesamples, nSpectra, nHarmonics);
  }

} //namespace astroaccelerate
