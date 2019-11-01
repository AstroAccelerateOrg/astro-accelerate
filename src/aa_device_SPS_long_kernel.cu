// Added by Karel Adamek 
#include "aa_device_SPS_long_kernel.hpp"

namespace astroaccelerate {

  //--------------------------------------------------------------------------------------
  //------------- device functions
  __device__ __inline__ float2 Calculate_SNR_const_memory(float x, float y, ushort Taps){
    float2 output;
    output.x = __fdividef( x , c_sqrt_taps[Taps] );
    output.y = __fdividef( y , c_sqrt_taps[Taps] );
    return(output);
  }

  __device__ __inline__ float2 Calculate_SNR(float x, float y, ushort Taps){
    float2 output;
    output.x = __frsqrt_rn(Taps)*x;
    output.y = __frsqrt_rn(Taps)*y;
    return(output);
  }

  __device__ __inline__ float2 Calculate_SNR_EACH(float x, float y, float StdDev){
    float2 output;
    output.x = __frsqrt_rn(StdDev)*x;
    output.y = __frsqrt_rn(StdDev)*y;
    return(output);
  }

  __device__ __inline__ void Compare( float2 *active_SNR, ushort2 *active_taps, float2 *passive_SNR, ushort *passive_taps){
    if(passive_SNR->x>active_SNR->x) {active_SNR->x=passive_SNR->x; active_taps->x=(*passive_taps);}
    if(passive_SNR->y>active_SNR->y) {active_SNR->y=passive_SNR->y; active_taps->y=(*passive_taps);}
  }

  __device__ __inline__ void Compare( float2 *active_SNR, ushort2 *active_taps, float2 *passive_SNR, ushort2 *passive_taps){
    if(passive_SNR->x>active_SNR->x) {active_SNR->x=passive_SNR->x; active_taps->x=passive_taps->x;}
    if(passive_SNR->y>active_SNR->y) {active_SNR->y=passive_SNR->y; active_taps->y=passive_taps->y;}
  }
  //------------- device functions
  //--------------------------------------------------------------------------------------


  __global__ void SPDT_GPU_1st_plane(float const* __restrict__ d_input, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float2 const* __restrict__ d_MSD, int nTimesamples, int nBoxcars, const int dtm) {
    __shared__ float2 s_input[PD_NTHREADS];
    __shared__ float2 s_SNRs[PD_NTHREADS];
    __shared__ ushort2 s_taps[PD_NTHREADS];
	
    int d, spos;
    size_t gpos;
    float Bw[8];
    float2 ftemp1, ftemp2, ftemp3, SNR;
    ushort2 taps, taps1, taps2;
    float2 stat;

	
    Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0; SNR.x=0; SNR.y=0;
	
    spos=blockIdx.x*(2*PD_NTHREADS - nBoxcars) + 2*threadIdx.x;
    gpos = (size_t)(blockIdx.y*nTimesamples) + (size_t)(spos);
	
    // Loading data and normalization
    if( (spos + 4)<nTimesamples ){
      ftemp1.x=d_input[gpos];	
      ftemp1.y=d_input[gpos + 1];
      ftemp2.x=d_input[gpos + 2];
      ftemp2.y=d_input[gpos + 3];
      ftemp3.x=d_input[gpos + 4];

      Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
      Bw[1]=ftemp1.x + ftemp1.y; Bw[5]=ftemp1.y + ftemp2.x;
      Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
      Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
    }
	
    spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
    if( threadIdx.x<(PD_NTHREADS - (nBoxcars>>1)) &&  spos<dtm) d_decimated[blockIdx.y*dtm + spos]=Bw[1];
		
    s_input[threadIdx.x].x = Bw[3];
    s_input[threadIdx.x].y = Bw[7];

    __syncthreads();
	
    stat = d_MSD[0];
    ftemp1.x = __fdividef(Bw[0] - stat.x, stat.y);
    ftemp1.y = __fdividef(Bw[4] - stat.x, stat.y);
    stat = d_MSD[1];
    ftemp2.x = __fdividef(Bw[1] - stat.x, stat.y);
    ftemp2.y = __fdividef(Bw[5] - stat.x, stat.y);
    if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=2;} else taps1.x=1;
    if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=2;} else taps1.y=1;
	
    stat = d_MSD[2];
    ftemp2.x = __fdividef(Bw[2] - stat.x, stat.y);
    ftemp2.y = __fdividef(Bw[6] - stat.x, stat.y);
    stat = d_MSD[3];
    ftemp3.x = __fdividef(Bw[3] - stat.x, stat.y);
    ftemp3.y = __fdividef(Bw[7] - stat.x, stat.y);
    if(ftemp3.x>ftemp2.x) {ftemp2.x=ftemp3.x; taps2.x=4;} else taps2.x=3;
    if(ftemp3.y>ftemp2.y) {ftemp2.y=ftemp3.y; taps2.y=4;} else taps2.y=3;
	
    if(ftemp1.x>ftemp2.x) {SNR.x=ftemp1.x; taps.x=taps1.x;} else {SNR.x=ftemp2.x; taps.x=taps2.x;}
    if(ftemp1.y>ftemp2.y) {SNR.y=ftemp1.y; taps.y=taps1.y;} else {SNR.y=ftemp2.y; taps.y=taps2.y;}

    s_SNRs[threadIdx.x]=SNR;
    s_taps[threadIdx.x]=taps;
	
    for(d=4; d<nBoxcars; d=d + 4){
      __syncthreads();
      spos=threadIdx.x-(d>>1);
      if(spos>=0){
	ftemp1 = s_input[spos];
	SNR = s_SNRs[spos];
	taps = s_taps[spos];
			
			
	Bw[0] = Bw[0] + ftemp1.x; Bw[4] = Bw[4] + ftemp1.y;
	Bw[1] = Bw[1] + ftemp1.x; Bw[5] = Bw[5] + ftemp1.y;
	Bw[2] = Bw[2] + ftemp1.x; Bw[6] = Bw[6] + ftemp1.y;
	Bw[3] = Bw[3] + ftemp1.x; Bw[7] = Bw[7] + ftemp1.y;
			
	stat = d_MSD[d];
	ftemp1.x = __fdividef( (Bw[0] - stat.x) , stat.y );
	ftemp1.y = __fdividef( (Bw[4] - stat.x) , stat.y );
	stat = d_MSD[d + 1];
	ftemp2.x = __fdividef( (Bw[1] - stat.x) , stat.y );
	ftemp2.y = __fdividef( (Bw[5] - stat.x) , stat.y );
	if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = d + 2;} else taps1.x = d + 1;
	if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = d + 2;} else taps1.y = d + 1;

	// Using constant memory
	stat = d_MSD[d + 2];
	ftemp2.x = __fdividef( (Bw[2] - stat.x) , stat.y );
	ftemp2.y = __fdividef( (Bw[6] - stat.x) , stat.y );
	stat = d_MSD[d + 3];
	ftemp3.x = __fdividef( (Bw[3] - stat.x) , stat.y );
	ftemp3.y = __fdividef( (Bw[7] - stat.x) , stat.y );
	if(ftemp3.x>ftemp2.x) {ftemp2.x = ftemp3.x; taps2.x = d + 4;} else taps2.x = d + 3;
	if(ftemp3.y>ftemp2.y) {ftemp2.y = ftemp3.y; taps2.y = d + 4;} else taps2.y = d + 3;
			
			
	if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = taps2.x;}
	if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = taps2.y;}
			
	if(ftemp1.x>SNR.x) {SNR.x = ftemp1.x; taps.x = taps1.x;}
	if(ftemp1.y>SNR.y) {SNR.y = ftemp1.y; taps.y = taps1.y;}
			
	s_SNRs[spos]=SNR;
	s_taps[spos]=taps;
      }
    }
	
    if(spos>=0) s_input[spos].x = Bw[3];
	
    __syncthreads();
	
    spos=blockIdx.x*(2*PD_NTHREADS - nBoxcars) + 2*threadIdx.x;
    if( threadIdx.x<(PD_NTHREADS - (nBoxcars>>1)) && spos<nTimesamples){
      gpos=(size_t)(blockIdx.y*nTimesamples) + (size_t)spos;
      d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
      d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
      d_output_taps[gpos] = s_taps[threadIdx.x].x;
      d_output_taps[gpos + 1] = s_taps[threadIdx.x].y;
		
      spos = blockIdx.x*(PD_NTHREADS - (nBoxcars>>1)) + threadIdx.x;
      if(spos<dtm) d_bv_out[blockIdx.y*dtm + spos] = s_input[threadIdx.x].x;
    }
  }


  __global__ void SPDT_GPU_Nth_plane(float const* __restrict__ d_input, float *d_bv_in, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float2 const* __restrict__ d_MSD, const int nTimesamples, const int nBoxcars, const int startTaps, const int DIT_value, const int dtm) {
    __shared__ float2 s_input[PD_NTHREADS];
    __shared__ float2 s_BV[PD_NTHREADS];
    __shared__ float2 s_SNRs[PD_NTHREADS];
    __shared__ ushort2 s_taps[PD_NTHREADS];
	
	
    int d, spos;
    size_t gpos;
    float Bw[8];
    float2 ftemp1, ftemp2, ftemp3, SNR;
    ushort2 taps, taps1, taps2;
    float2 stat;
	
    Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0;
    s_BV[threadIdx.x].x=0; s_BV[threadIdx.x].y=0;
	
    spos=blockIdx.x*(2*PD_NTHREADS - nBoxcars) + 2*threadIdx.x;
    gpos=(size_t)(blockIdx.y*nTimesamples) + (size_t)spos;
	
    if( (spos + 4)<nTimesamples ){
      ftemp1.x=d_input[gpos];
      ftemp1.y=d_input[gpos + 1];
      ftemp2.x=d_input[gpos + 2];
      ftemp2.y=d_input[gpos + 3];
      ftemp3.x=d_input[gpos + 4];
		
		
      Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
      Bw[1]=ftemp1.x + ftemp1.y; Bw[5]=ftemp1.y + ftemp2.x;
      Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
      Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
		
      s_BV[threadIdx.x].x=d_bv_in[gpos];//   - BV_mean;	
      s_BV[threadIdx.x].y=d_bv_in[gpos + 1];// - BV_mean;
    }
	
    spos = blockIdx.x*(PD_NTHREADS - (nBoxcars>>1)) + threadIdx.x;
    if( threadIdx.x<(PD_NTHREADS - (nBoxcars>>1)) &&  spos<dtm) d_decimated[blockIdx.y*dtm + spos]=Bw[1];
	
    s_input[threadIdx.x].x = Bw[3];
    s_input[threadIdx.x].y = Bw[7];
	
    SNR.x=0; SNR.y=0;
    __syncthreads();
	
    taps1.x = 1; taps1.y=2;
    stat = d_MSD[0];
    ftemp1.x = __fdividef( (s_BV[threadIdx.x].x + (Bw[0] - stat.x)) , stat.y );
    ftemp1.y = __fdividef( (s_BV[threadIdx.x].y + (Bw[4] - stat.x)) , stat.y );
    stat = d_MSD[1];
    ftemp2.x = __fdividef( (s_BV[threadIdx.x].x + (Bw[1] - stat.x)) , stat.y );
    ftemp2.y = __fdividef( (s_BV[threadIdx.x].y + (Bw[5] - stat.x)) , stat.y );
    if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=2;} else taps1.x=1;
    if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=2;} else taps1.y=1;

    taps2.x=3; taps2.y=4;
    stat = d_MSD[2];
    ftemp2.x = __fdividef( (s_BV[threadIdx.x].x + (Bw[2] - stat.x)) , stat.y );
    ftemp2.y = __fdividef( (s_BV[threadIdx.x].y + (Bw[6] - stat.x)) , stat.y );
    stat = d_MSD[3];
    ftemp3.x = __fdividef( (s_BV[threadIdx.x].x + (Bw[3] - stat.x)) , stat.y );
    ftemp3.y = __fdividef( (s_BV[threadIdx.x].y + (Bw[7] - stat.x)) , stat.y );
    if(ftemp3.x>ftemp2.x) {ftemp2.x=ftemp3.x; taps2.x=4;} else taps2.x=3;
    if(ftemp3.y>ftemp2.y) {ftemp2.y=ftemp3.y; taps2.y=4;} else taps2.y=3;
	
    if(ftemp1.x>ftemp2.x) {SNR.x=ftemp1.x; taps.x=taps1.x;} else {SNR.x=ftemp2.x; taps.x=taps2.x;}
    if(ftemp1.y>ftemp2.y) {SNR.y=ftemp1.y; taps.y=taps1.y;} else {SNR.y=ftemp2.y; taps.y=taps2.y;}
	
    s_SNRs[threadIdx.x]=SNR;
    s_taps[threadIdx.x]=taps;
	
    for(d=4; d<nBoxcars; d=d + 4){
      __syncthreads();
      spos=threadIdx.x-(d>>1);
      if(spos>=0){
	ftemp1 = s_input[spos];
	SNR = s_SNRs[spos];
	taps = s_taps[spos];
			
	Bw[0] = Bw[0] + ftemp1.x; Bw[4] = Bw[4] + ftemp1.y;
	Bw[1] = Bw[1] + ftemp1.x; Bw[5] = Bw[5] + ftemp1.y;
	Bw[2] = Bw[2] + ftemp1.x; Bw[6] = Bw[6] + ftemp1.y;
	Bw[3] = Bw[3] + ftemp1.x; Bw[7] = Bw[7] + ftemp1.y;
			
	taps1.x = d + 1; taps1.y=d + 2;
	stat = d_MSD[d];
	ftemp1.x = __fdividef( (s_BV[spos].x + (Bw[0] - stat.x)) , stat.y );
	ftemp1.y = __fdividef( (s_BV[spos].y + (Bw[4] - stat.x)) , stat.y );
	stat = d_MSD[d + 1];
	ftemp2.x = __fdividef( (s_BV[spos].x + (Bw[1] - stat.x)) , stat.y );
	ftemp2.y = __fdividef( (s_BV[spos].y + (Bw[5] - stat.x)) , stat.y );
	if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = d + 2;} else taps1.x = d + 1;
	if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = d + 2;} else taps1.y = d + 1;
			
	taps2.x=d + 3; taps2.y=d + 4;
	stat = d_MSD[d + 2];
	ftemp2.x = __fdividef( (s_BV[spos].x + (Bw[2] - stat.x)) , stat.y);
	ftemp2.y = __fdividef( (s_BV[spos].y + (Bw[6] - stat.x)) , stat.y);
	stat = d_MSD[d + 3];
	ftemp3.x = __fdividef( (s_BV[spos].x + (Bw[3] - stat.x)) , stat.y);
	ftemp3.y = __fdividef( (s_BV[spos].y + (Bw[7] - stat.x)) , stat.y);
	if(ftemp3.x>ftemp2.x) {ftemp2.x = ftemp3.x; taps2.x = d + 4;} else taps2.x = d + 3;
	if(ftemp3.y>ftemp2.y) {ftemp2.y = ftemp3.y; taps2.y = d + 4;} else taps2.y = d + 3;
			
	if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = taps2.x;}
	if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = taps2.y;}
			
	if(ftemp1.x>SNR.x) {SNR.x = ftemp1.x; taps.x = taps1.x;}
	if(ftemp1.y>SNR.y) {SNR.y = ftemp1.y; taps.y = taps1.y;}
			
	s_SNRs[spos]=SNR;
	s_taps[spos]=taps;
      }

    }
	
    if(spos>=0) s_input[spos].x = s_BV[spos].x + Bw[3];
	
    __syncthreads();
	
    spos = blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
    if( threadIdx.x<(PD_NTHREADS - (nBoxcars>>1))  && spos<(nTimesamples) ){
      gpos=(size_t)(blockIdx.y*nTimesamples) + (size_t)spos;
      d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
      d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
      d_output_taps[gpos] = startTaps + s_taps[threadIdx.x].x*DIT_value;
      d_output_taps[gpos + 1] = startTaps + s_taps[threadIdx.x].y*DIT_value;
		
      spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
      if(spos<dtm) d_bv_out[blockIdx.y*dtm + spos] = s_input[threadIdx.x].x;
    }
  }

  /** \brief Kernel wrapper function to SPDT_GPU_1st_plane kernel function. */
  void call_kernel_SPDT_GPU_1st_plane(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_bv_out,
				      float *const d_decimated,
				      float *const d_output_SNR, ushort *const d_output_taps, float2 const *const d_MSD,
				      const int &nTimesamples, const int &nBoxcars, const int &dtm) {
    SPDT_GPU_1st_plane<<<grid_size,block_size>>>(d_input, d_bv_out,
						 d_decimated,
						 d_output_SNR, d_output_taps, d_MSD,
						 nTimesamples, nBoxcars, dtm);
  }

  /** \brief Kernel wrapper function to SPDT_GPU_Nth_plane kernel function. */
  void call_kernel_SPDT_GPU_Nth_plane(const dim3 &grid_size, const dim3 &block_size,
				      float const *const d_input, float *const d_bv_in, float *const d_bv_out,
				      float *const d_decimated, float *const d_output_SNR, ushort *const d_output_taps,
				      float2 const *const d_MSD, const int &nTimesamples,
				      const int &nBoxcars, const int &startTaps, const int &DIT_value, const int &dtm) {
    SPDT_GPU_Nth_plane<<<grid_size,block_size>>>(d_input, d_bv_in, d_bv_out,
						 d_decimated, d_output_SNR, d_output_taps,
						 d_MSD, nTimesamples,
						 nBoxcars, startTaps, DIT_value, dtm);
  }

} //namespace astroaccelerate
