#include "device_zero_dm_outliers_kernel.hpp"
#include "params.hpp"
#include <stdio.h>
#include <time.h>

//{{{ zero_dm

void zero_dm_outliers(unsigned short* d_input, int nchans, int nsamp) {

  int divisions_in_t = 224;
  int num_blocks_t   = nsamp / divisions_in_t;

  printf("\nZDM OUTLIERS!");
  printf("\n%d %d", nsamp, nchans);
  printf("\n%d %d", divisions_in_t, 1);
  printf("\n%d %d", num_blocks_t, 1);

  dim3 threads_per_block(divisions_in_t, 1);
  dim3 num_blocks(num_blocks_t, 1);

  clock_t start_t, end_t;
  start_t = clock();

  call_kernel_zero_dm_outliers_kernel_one(
      num_blocks, threads_per_block, d_input, nchans, nsamp);
  cudaDeviceSynchronize();

  int divisions_in_c = 100;
  int num_blocks_c   = nchans / divisions_in_c;

  printf("\nZDM OUTLIERS!");
  printf("\n%d %d", nsamp, nchans);
  printf("\n%d %d", divisions_in_c, 1);
  printf("\n%d %d", num_blocks_c, 1);

  dim3 threads_per_block_c(divisions_in_c, 1);
  dim3 c_blocks(num_blocks_c, 1);

  call_kernel_zero_dm_outliers_kernel_two(
      c_blocks, threads_per_block_c, d_input, nchans, nsamp);

  end_t       = clock();
  double time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("\nPerformed ZDM: %lf (GPU estimate)", time);

  //}}}
}

//}}}
