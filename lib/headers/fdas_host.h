/* FDAS host functions header */
#ifndef ASTROACCELERATE_FDAS_HOST_H_
#define ASTROACCELERATE_FDAS_HOST_H_

#include <cufft.h>
#include "fdas_device.h"
#include "fdas_util.h"
#include "params.h"

#define LOWACC 0
#define SLIGHT 299792458.0

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

typedef struct{
  cufftHandle realplan;
  cufftHandle forwardplan;
  cufftHandle inverseplan;
}fdas_cufftplan;

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

//external presto functions
//------------------------------
extern "C" {
  cufftComplex *presto_gen_z_response( double z, int numkern, int numbetween);

  void presto_place_complex_kernel(cufftComplex * kernel, int numkernel, cufftComplex * result, int numresult);

  int presto_z_resp_halfwidth(double z, int accuracy);

  void presto_dered_sig(cufftComplex * fft,  unsigned long numamps);

  void presto_norm(cufftComplex * fft,  unsigned long numamps);

  double candidate_sigma(double power, int numsum, double numtrials);
}

#endif /* FDAS_HOST_H */
