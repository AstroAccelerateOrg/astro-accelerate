#include <cuda_runtime.h>

#include "aa_bin_gpu.hpp"
#include "aa_gpu_timer.hpp"

namespace astroaccelerate {

void bin_gpu(unsigned short *const d_input, float *const d_output, const int nchans, const int nsamp) {
    
    int divisions_in_t = BINDIVINT;
    int divisions_in_f = BINDIVINF;
    int num_blocks_t = (int) ( ( nsamp + 1 ) / ( 2 * divisions_in_t ) );
    int num_blocks_f = nchans / divisions_in_f;
    
    dim3 threads_per_block(divisions_in_t, divisions_in_f);
    dim3 num_blocks(num_blocks_t, num_blocks_f);
    //printf("\nDIVISOR:\t%f", (float)(nsamp)/(2*divisions_in_t));
    //printf("\ndivisions_in_t:%d\tdivisions_in_f:%d",divisions_in_t, divisions_in_f);
    //printf("\nnum_blocks_t:%d\tnum_blocks_f:%d\n",num_blocks_t,num_blocks_f);

    size_t size_chunk = (size_t)(nchans)*(size_t)(nsamp);
    cudaMemset(d_output, 0, size_chunk*sizeof(float));
    call_kernel_bin(num_blocks, threads_per_block, d_input, d_output, nsamp);
    
    int swap_divisions_in_t = CT;
    int swap_divisions_in_f = CF;
    int swap_num_blocks_t = nsamp / swap_divisions_in_t;
    int swap_num_blocks_f = nchans / swap_divisions_in_f;
    
    dim3 swap_threads_per_block(swap_divisions_in_t, swap_divisions_in_f);
    dim3 swap_num_blocks(swap_num_blocks_t, swap_num_blocks_f);

    cudaDeviceSynchronize();
    call_kernel_swap(swap_num_blocks, swap_threads_per_block, d_input, d_output, nchans, nsamp);
    cudaDeviceSynchronize();
    
    // end_t = clock();
    // double time = double ( end_t - start_t ) / CLOCKS_PER_SEC;
    //printf("\nPerformed Bin: %f (GPU estimate)", time);
    //printf("\nGops based on %.2f ops per channel per tsamp: %f",14.0,((15.0*(divisions_in_t*divisions_in_f*num_blocks_t*num_blocks_f))/(time))/1000000000.0);
    //printf("\nBN Device memory bandwidth in GB/s: %f", (2*(sizeof(float)+sizeof(unsigned short))*(divisions_in_t*divisions_in_f*num_blocks_t*num_blocks_f))/(time)/1000000000.0);
//    size_chunk = (size_t)(nchans)*(size_t)(nsamp);
    cudaMemset(d_output, 0, size_chunk*sizeof(float));
}

int GPU_DiT_v2_wrapper(float *d_input, float *d_output, int nDMs, int nTimesamples) {
  aa_gpu_timer timer;

  //---------> Task specific
  int ut;

  //---------> CUDA block and CUDA grid parameters
  int nThreads = 2*WARP;
  int nCUDAblocks_x=nTimesamples/(DIT_ELEMENTS_PER_THREAD*nThreads*2);
  int nRest = nTimesamples%(DIT_ELEMENTS_PER_THREAD*nThreads*2);
  if(nRest>0) nCUDAblocks_x++;
  int nCUDAblocks_y=nDMs/DIT_YSTEP;

  dim3 gridSize(nCUDAblocks_x, nCUDAblocks_y, 1);
  dim3 blockSize(nThreads, DIT_YSTEP, 1);

  if( (nTimesamples>>1)==0 ) return(-1);

  // ----------------------------------------------->
  // --------> Measured part (Pulse detection FIR)
  //---------> Pulse detection FIR
  call_kernel_DiT_GPU_v2(gridSize,
			 blockSize,
			 d_input,
			 d_output,
			 nDMs,
			 nTimesamples,
			 (nTimesamples>>1));
  // --------> Measured part (Pulse detection FIR)
  // ----------------------------------------------->

  ut=0;
  return(ut);
}


} //namespace astroaccelerate
