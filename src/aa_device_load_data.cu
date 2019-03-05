#include "aa_device_load_data.hpp"
#include "aa_device_dedispersion_kernel.hpp"
#include "aa_device_SPS_inplace_kernel.hpp"

#include "aa_params.hpp"

namespace astroaccelerate {

  /**
   * \brief Function to load data from host memory into GPU memory.
   * \warning If the file extension of this file is *.cpp, then the code will compile but there will be a runtime CUDA error when copying to device memory.
   */
  
  void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short const*const host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts) {
    if(i == -1) {
      const long int length = ( t_processed + maxshift );
      const size_t size = (size_t)nchans * (size_t)length * (size_t)sizeof(unsigned short);
      //checkCudaErrors(cudaGetLastError());
      cudaMemcpy(device_pointer, host_pointer, size, cudaMemcpyHostToDevice);
      //checkCudaErrors(cudaGetLastError());
      set_device_constants_dedispersion_kernel(nchans, length, t_processed, dmshifts);
    }
    else if(i > 0) {
      const long int length = ( t_processed + maxshift );
      set_device_constants_dedispersion_kernel(length, t_processed);
    }

    float h_sqrt_taps[PD_MAXTAPS + 1];
    for(int f = 0; f <= PD_MAXTAPS; f++) {
      h_sqrt_taps[f] = (float) sqrt((double) f);
    }
    cudaMemcpyToSymbol(c_sqrt_taps, &h_sqrt_taps, ( PD_MAXTAPS + 1 ) * sizeof(float));
  }

} //namespace astroaccelerate
