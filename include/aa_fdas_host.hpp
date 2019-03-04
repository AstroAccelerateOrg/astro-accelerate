#ifndef ASTRO_ACCELERATE_AA_FDAS_HOST_HPP
#define ASTRO_ACCELERATE_AA_FDAS_HOST_HPP

#include <cufft.h>
#include "aa_fdas_device.hpp"
#include "aa_fdas_util.hpp"
#include "aa_params.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include "presto_funcs.hpp"
#include "presto.hpp"
#include "aa_params.hpp"
#include "aa_fdas_test_parameters.hpp"
//#include <helper_functions.h>
#include <curand.h>
#include <libgen.h>
//#include <random> // C++11 to use normal distribution

namespace astroaccelerate {

#define LOWACC 0
#define SLIGHT 299792458.0

  /**
   * \struct fdas_new_acc_sig
   * \brief Struct to hold new signal data for fdas.
   */
  typedef struct{
    float* acc_signal;
    int nsamps;
    double freq0;
    int mul;
    int zval; 
    int nsig;
    int nharms;
    double duty;
    float sigamp;
  }fdas_new_acc_sig;

  /**
   * \struct fdas_gpuarrays
   * \brief Struct to hold the device data.
   */
  typedef struct{
    float* d_in_signal; 
    float2* d_fft_signal;
    float2  *d_ext_data;
    float2 *d_kernel;
    float *d_ffdot_pwr;
    float *d_ffdot_summed;
    float2 *d_ffdot_cpx;
    float2 *ip_edge_points;// edge points for interbinning in kfft
    float *d_fdas_peak_list; // added by KA
    size_t mem_insig;
    size_t mem_rfft;
    size_t mem_extsig;
    size_t mem_ffdot;
    size_t mem_ffdot_cpx;
    size_t mem_ipedge; // edge points for interbinning in kfft
    size_t mem_max_list_size; // maximum length of candate list in bytes added by KA
  }fdas_gpuarrays;

  /**
   * \struct fdas_cufftplan
   * \brief Type to hold cufft plan data.
   */
  typedef struct{
    cufftHandle realplan;
    cufftHandle forwardplan;
    cufftHandle inverseplan;
  }fdas_cufftplan;

  /**
   * \struct fdas_params
   * \brief Struct to hold fdas parameter metadata.
   */
  typedef struct{
    int  nsamps;
    int rfftlen;
    int sigblock;
    int nblocks; 
    int offset; 
    int siglen;
    int extlen;
    int max_list_length; // maximum number of rows in the list
    unsigned int ffdotlen;
    unsigned int ffdotlen_cpx;
    float scale;
    float tsamp;
  }fdas_params;

  //function declarations

  void fdas_print_params_h();

  void fdas_cuda_check_devices(int devid);

  void fdas_alloc_gpu_arrays(fdas_gpuarrays *arrays,  cmd_args *cmdargs);

  void fdas_free_gpu_arrays(fdas_gpuarrays *arrays,  cmd_args *cmdargs);

  void fdas_create_acc_sig(fdas_new_acc_sig *acc_sig, cmd_args *cmdargs);

  void fdas_create_acc_kernels(cufftComplex* d_kernel, cmd_args *cmdargs );

  void fdas_cuda_create_fftplans( fdas_cufftplan *fftplans, fdas_params *params);

  void fdas_cuda_basic(fdas_cufftplan *fftplans, fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params );

  void fdas_cuda_customfft(fdas_cufftplan *fftplans, fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params );

  void fdas_write_list(fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params, float *h_MSD, float dm_low, int dm_count, float dm_step, unsigned int list_size);

  void fdas_write_ffdot(fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params, float dm_low, int dm_count, float dm_step );

  void fdas_write_test_ffdot(fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params, float dm_low, int dm_count, float dm_step );

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_FDAS_HOST_HPP
