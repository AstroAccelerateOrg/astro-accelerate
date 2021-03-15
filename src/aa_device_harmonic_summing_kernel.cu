#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_params.hpp"
#include "aa_device_harmonic_summing_kernel.hpp"

namespace astroaccelerate {

__global__ void simple_harmonic_sum_GPU_kernel(float const* __restrict__ d_input, float *d_output_SNR, ushort *d_output_harmonics, float *d_MSD, int nTimesamples, int nSpectra, int nHarmonics){
  float HS_value, temp_SNR, SNR;
  ushort max_SNR_harmonic;
  int pos;

  // reading 0th harmonic, i.e. fundamental frequency
  pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
  if( (blockIdx.y*blockDim.x + threadIdx.x)<nSpectra ){
    HS_value = __ldg(&d_input[pos]);
    SNR = (HS_value - __ldg(&d_MSD[0]))/(__ldg(&d_MSD[1]));
    max_SNR_harmonic = 0;
    
    if(blockIdx.x>0) {
      for(int f=1; f<nHarmonics; f++) {
        if( (blockIdx.x + f*blockIdx.x)<nTimesamples ) {
          pos = (blockIdx.x + f*blockIdx.x)*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
          HS_value = HS_value + __ldg(&d_input[pos]);
          temp_SNR = (HS_value - __ldg(&d_MSD[f*2]))/(__ldg(&d_MSD[2*f+1])); //assuming white noise 
          if(temp_SNR > SNR) {
            SNR = temp_SNR;
            max_SNR_harmonic = f;
          }
        }
      }
    }
    
    pos = blockIdx.x*nSpectra + blockIdx.y*blockDim.x + threadIdx.x;
    d_output_SNR[pos] = SNR;
    d_output_harmonics[pos] = max_SNR_harmonic;
  }
}

__inline__ __device__ float remove_scalloping_loss(float Xm2, float Xm1, float X0, float Xp1, float Xp2){
    return(X0 + (1.88494/2.0)*(Xm1 + Xp1) + (0.88494/2.0)*(Xm2 + Xp2));
}

template<class const_params>
__inline__ __device__ void get_frequency_bins(float *down, float *step, float const* __restrict__ d_input, int pos){
    if(const_params::remove_scalloping_loss) {
        float Xm2 = d_input[pos - 2];
        float Xm1 = d_input[pos - 1];
        float X0  = d_input[pos + 0];
        float Xp1 = d_input[pos + 1];
        float Xp2 = d_input[pos + 2];
        float Xp3 = d_input[pos + 3];
        (*down) = remove_scalloping_loss(Xm2, Xm1, X0, Xp1, Xp2);
        (*step) = remove_scalloping_loss(Xm1, X0, Xp1, Xp2, Xp3);
    }
    else {
        (*down) = d_input[pos];
        (*step) = d_input[pos + 1];
    }
}

template<class const_params>
__global__ void greedy_harmonic_sum_GPU_kernel(float *d_maxSNR, ushort *d_maxHarmonics, float const* __restrict__ d_input, float const* __restrict__ d_MSD, int nTimesamples, int nDMs, int nHarmonics){
    __shared__ float s_MSD[64];
    float SNR;
    float partial_sum, maxSNR;
    int maxHarmonics;
    
    if(threadIdx.x<nHarmonics) {
        s_MSD[threadIdx.x]      = d_MSD[2*threadIdx.x];
        s_MSD[32 + threadIdx.x] = d_MSD[2*threadIdx.x + 1];
    }
    
    __syncthreads();
    
    int data_shift=0;
    float data_down = 0, data_step = 0;
    int pos = GHRMS_NTHREADS*blockIdx.x + threadIdx.x;
    if( pos > 1 && pos + 3 < nTimesamples) {
        int block_pos = blockIdx.y*nTimesamples + pos;
        get_frequency_bins<const_params>(&data_down, &data_step, d_input, block_pos);
    }
    if( data_step > data_down ) {
        data_shift++;
        partial_sum = data_step;
    }
    else {
        partial_sum = data_down;
    }
    maxSNR = fdividef( (partial_sum - s_MSD[0]), s_MSD[1]);
    maxHarmonics = 0;
    
    for(int h=1; h<nHarmonics; h++){
        int pos = (h+1)*(GHRMS_NTHREADS*blockIdx.x + threadIdx.x) + data_shift;
        float data_down = 0, data_step = 0;
        if( pos > 1 && pos + 3 < nTimesamples) {
            int block_pos = blockIdx.y*nTimesamples + pos;
            get_frequency_bins<const_params>(&data_down, &data_step, d_input, block_pos);
        }
        
        if( data_step > data_down ) {
            data_shift++;
            partial_sum = partial_sum + data_step;
        }
        else {
            partial_sum = partial_sum + data_down;
        }
        
        SNR = fdividef( (partial_sum - s_MSD[2*h]), s_MSD[2*h+1] );
        if( SNR > maxSNR ) {
            maxSNR = SNR;
            maxHarmonics = h;
        }
    }
    
    pos = GHRMS_NTHREADS*blockIdx.x + threadIdx.x;
    if(pos < nTimesamples) {
        d_maxSNR[blockIdx.y*nTimesamples + pos] = maxSNR;
        d_maxHarmonics[blockIdx.y*nTimesamples + pos] = (ushort) maxHarmonics;
    }
}


template<class const_params>
__inline__ __device__ void get_frequency_bin_value(float *frequency_bin, float const* __restrict__ data, int pos){
    if(const_params::remove_scalloping_loss){
        float Xm2 = data[pos - 2];
        float Xm1 = data[pos - 1];
        float X0  = data[pos];
        float Xp1 = data[pos + 1];
        float Xp2 = data[pos + 2];
        (*frequency_bin) = X0 + (1.88494/2.0)*(Xm1 + Xp1) + (0.88494/2.0)*(Xm2 + Xp2);
    }
    else {
        (*frequency_bin) = data[pos];
    }
}


template<class const_params>
__global__ void presto_harmonic_sum_GPU_kernel(float *d_maxSNR, ushort *d_maxHarmonics, float const* __restrict__ d_input, float const* __restrict__ d_MSD, int nTimesamples, int nDMs, int nHarmonics){
    __shared__ float s_MSD[64];
    float SNR;
    float partial_sum, maxSNR, frequency_bin, fundamental;
    int maxHarmonics, pos;
    
    if(threadIdx.x<nHarmonics) {
        s_MSD[2*threadIdx.x]   = d_MSD[2*threadIdx.x];
        s_MSD[2*threadIdx.x+1] = d_MSD[2*threadIdx.x+1];
    }
    
    __syncthreads();
    
    partial_sum = 0;
    frequency_bin = 0;
    fundamental = 0;
    pos = GHRMS_NTHREADS*blockIdx.x + threadIdx.x;
    if( (pos > 1) && (pos + 2) < nTimesamples ) {
        int block_pos = blockIdx.y*nTimesamples + pos;
        get_frequency_bin_value<const_params>(&fundamental, d_input, block_pos);
    }
    partial_sum = fundamental;
    maxSNR = fdividef( (partial_sum - s_MSD[0]), s_MSD[1]);
    maxHarmonics = 0;
    
    if( pos > 1 && (pos + 2) < nTimesamples ) {
        for(int i = 1; i < nHarmonics; i++){ //i + 1 = num. of harmonic added;
            partial_sum = fundamental;
            double fundamental_fraction = ((double) pos)/((double) (i + 1));
            for(int f= 1; f<=i; f++){
                int new_pos = (int) ( ( ((double) f)*fundamental_fraction ) + 0.5 );
                int block_pos = blockIdx.y*nTimesamples + new_pos;
                frequency_bin = 0;
                if( new_pos > 1 && (new_pos + 2) < nTimesamples ) {
                    get_frequency_bin_value<const_params>(&frequency_bin, d_input, block_pos);
                }
                partial_sum = partial_sum + frequency_bin;
            }
            SNR = fdividef( (partial_sum - s_MSD[2*i]), s_MSD[2*i + 1]);
            if(SNR>maxSNR) {
                maxSNR = SNR;
                maxHarmonics = i-1;
            }
        }
    }
    //----------------------------------------------
    
    __syncthreads();
    
    pos = GHRMS_NTHREADS*blockIdx.x + threadIdx.x;
    if( pos < nTimesamples ){
        int block_pos = blockIdx.y*nTimesamples + pos;
        d_maxSNR[block_pos] = maxSNR;
        d_maxHarmonics[block_pos] = maxHarmonics;
    }
}





//-------------------------------------------------------------------
//------------------------------- Callers ---------------------------
//-------------------------------------------------------------------

  /** \brief Kernel wrapper function for simple_harmonic_sum_GPU_kernel kernel function. */
  void call_kernel_simple_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nSpectra,
      const int &nHarmonics
  ) {
    simple_harmonic_sum_GPU_kernel<<<grid_size, block_size>>>(d_input, d_output_SNR, d_output_harmonics, d_MSD, nTimesamples, nSpectra, nHarmonics);
  }

  /** \brief Kernel wrapper function for call_kernel_greedy_harmonic_sum_GPU_kernel kernel function. */
  void call_kernel_greedy_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nDMs,
      const int &nHarmonics,
      bool enable_scalloping_loss_removal
  ) {
    if(enable_scalloping_loss_removal) {
      greedy_harmonic_sum_GPU_kernel<HRMS_remove_scalloping_loss><<<grid_size, block_size>>>(
          d_output_SNR,
          d_output_harmonics,
          d_input,
          d_MSD,
          nTimesamples,
          nDMs,
          nHarmonics
      );
    }
    else {
      greedy_harmonic_sum_GPU_kernel<HRMS_normal><<<grid_size, block_size>>>(
          d_output_SNR,
          d_output_harmonics,
          d_input,
          d_MSD,
          nTimesamples,
          nDMs,
          nHarmonics
      );
    }
  }

  /** \brief Kernel wrapper function for presto_harmonic_sum_GPU_kernel kernel function. */
  void call_kernel_presto_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nDMs,
      const int &nHarmonics,
      bool enable_scalloping_loss_removal
  ) {
    if(enable_scalloping_loss_removal) {
      presto_harmonic_sum_GPU_kernel<HRMS_remove_scalloping_loss><<<grid_size, block_size>>>(
          d_output_SNR,
          d_output_harmonics,
          d_input,
          d_MSD,
          nTimesamples,
          nDMs,
          nHarmonics
      );
    }
    else {
      presto_harmonic_sum_GPU_kernel<HRMS_normal><<<grid_size, block_size>>>(
          d_output_SNR,
          d_output_harmonics,
          d_input,
          d_MSD,
          nTimesamples,
          nDMs,
          nHarmonics
      );
    }
  }
} //namespace astroaccelerate




