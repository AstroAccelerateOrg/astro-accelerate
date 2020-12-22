#include "aa_zero_dm.hpp"

namespace astroaccelerate {

  /** \brief Performs zero_dm. */
  void zero_dm(unsigned short *const d_input, const int nchans, const int nsamp, const int nbits) {
    
    int threads_for_sum  = 256;
    int num_blocks_t    = nsamp;
    int shared_memory = threads_for_sum*4;
    
    //printf("\nCORNER TURN!");
    //printf("\n%d %d %d", nsamp, nchans, threads_for_sum);
    //printf("\n%d %d", CT, 1);
    //printf("\n%d %d", num_blocks_t, 1);
    
    dim3 threads_per_block(threads_for_sum, 1);
    dim3 num_blocks(num_blocks_t,1);
    
    clock_t start_t, end_t;
    start_t = clock();
    
    float normalization_factor = ((pow(2,nbits)-1)/2);
    call_kernel_zero_dm_kernel(num_blocks, threads_per_block, shared_memory, d_input, nchans, nsamp, nbits, normalization_factor);
    cudaDeviceSynchronize();
    
    end_t = clock();
    double time = (double)(end_t-start_t) / CLOCKS_PER_SEC;
    printf("\nPerformed ZDM: %lf (GPU estimate)", time);    
  }
} //namespace astroaccelerate
