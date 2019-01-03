#ifndef ASTRO_ACCELERATE_AA_DEVICE_MSD_CONFIGURATION_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_MSD_CONFIGURATION_HPP

#include "aa_params.hpp"
#include <stdio.h>
#include <cuda_runtime.h>

namespace astroaccelerate {

  /**
   * \class MSD_Configuration aa_device_MSD_Configuration.hpp "include/aa_device_MSD_Configuration.hpp"
   * \brief Class that AstroAccelerate uses to manage the Mean Standard Deviation (MSD) configuration.
   * \brief Users should not use this class directly.
   * \author -
   * \date -
   */
class MSD_Configuration {
public:
  ushort2 nBlocks;
  ushort2 nSteps;
  dim3 partials_gridSize;
  dim3 partials_blockSize;
  dim3 final_gridSize;
  dim3 final_blockSize;
  int nBlocks_total;
  int address;
  size_t nTimesamples;
  size_t nDMs;
  size_t offset;
	

  int Choose_Divider(int number, int max_divider) {
    int seive[12]={2, 3, 4, 5, 7, 11, 13, 17, 19, 23, 29, 31};
    int f, nRest, nBlocks, N, N_accepted;

    N=1; N_accepted=1;
    do {
      N=1;
      for(f=0; f<12; f++) {
	nBlocks=number/seive[f];
	nRest=number - nBlocks*seive[f];
	if(nRest==0) {
	  N=seive[f];
	  N_accepted=N_accepted*N;
	  break;
	}
      }
      number=number/N;
    } while((N_accepted)<=max_divider && N>1);

    return(N_accepted/N);
  }

  int ChooseNumberOfThreads() {
    int nThreads=2048;
    int itemp=0;

    while(itemp==0 && nThreads>32) {
      nThreads=(nThreads>>1);
      itemp=(int)(nBlocks_total/(nThreads*32));
    }
    if(nThreads<32) nThreads=32;

    return(nThreads);
  }

  void Calculate_Kernel_Parameters() {
    int nRest;
    int nThreads;

    nSteps.x  = PD_NTHREADS;
    nBlocks.x = (int)((nTimesamples-offset)/(nSteps.x));
    nRest     = (nTimesamples - offset) - nBlocks.x*nSteps.x;
    if(nRest>0) nBlocks.x++;

    nSteps.y  = Choose_Divider(nDMs, 64);
    nBlocks.y = nDMs/nSteps.y;

    nBlocks_total = nBlocks.x*nBlocks.y;
    nThreads = ChooseNumberOfThreads();

    partials_gridSize.x = nBlocks.x;
    partials_gridSize.y = nBlocks.y;
    partials_gridSize.z = 1;

    partials_blockSize.x = PD_NTHREADS;
    partials_blockSize.y = 1;
    partials_blockSize.z = 1;

    final_gridSize.x = 1;
    final_gridSize.y = 1;
    final_gridSize.z = 1;

    final_blockSize.x = nThreads;
    final_blockSize.y = 1;
    final_blockSize.z = 1;
  }

  void reset() {
    nBlocks.x = 0; nBlocks.y = 0;
    nSteps.y = 0; nSteps.y = 0;
    partials_gridSize.x = 1; partials_gridSize.y = 1; partials_gridSize.z = 1;
    partials_blockSize.x = 1; partials_blockSize.y = 1; partials_blockSize.z = 1;
    final_gridSize.x = 1; final_gridSize.y = 1; final_gridSize.z = 1;
    final_blockSize.x = 1; final_blockSize.y = 1; final_blockSize.z = 1;
    nBlocks_total = 0;
    address = 0;
    nTimesamples = 0;
    nDMs = 0;
    offset = 0;
  }

  MSD_Configuration(void) {
    reset();
  }

  MSD_Configuration(size_t t_nTimesamples, size_t t_nDMs, size_t t_offset, int t_address) {
    nTimesamples = t_nTimesamples;
    nDMs         = t_nDMs;
    offset       = t_offset;
    address      = t_address;
		
    Calculate_Kernel_Parameters();
  }

  void print() {
    printf("         nTimesamples:%zu; nDMs:%zu; offset:%zu;\n", nTimesamples, nDMs, offset);
    printf("         nBlocks:[%d;%d]; nSteps:[%d;%d]; nBlocks_total:%d; address:%d;\n", nBlocks.x, nBlocks.y, nSteps.x, nSteps.y, nBlocks_total, address);
    printf("         partials_gridSize=[%d;%d;%d]; partials_blockSize=[%d;%d;%d]\n", partials_gridSize.x, partials_gridSize.y, partials_gridSize.z, partials_blockSize.x, partials_blockSize.y, partials_blockSize.z);
    printf("         final_gridSize=[%d;%d;%d]; final_blockSize=[%d;%d;%d]\n", final_gridSize.x, final_gridSize.y, final_gridSize.z, final_blockSize.x, final_blockSize.y, final_blockSize.z);
  }
};

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_DEVICE_MSD_CONFIGURATION_HPP
