// Added by Karel Adamek 

#ifndef SPS_LONG_KERNEL_H_
#define SPS_LONG_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

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




//******************************************************************************************
//******************************************************************************************
//******************************************************************************************

__global__ void PD_GPU_1st_BLN(float const* __restrict__ d_input, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int nTimesamples, int nBoxcars, const int dtm) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, ftemp4, SNR;
	ushort2 taps, taps1, taps3;
	ushort Taps;
	float signal_mean=d_MSD[0];
	float signal_sd=d_MSD[1];
	
	Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0; SNR.x=0; SNR.y=0;
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	gpos = blockIdx.y*nTimesamples + spos;
	
	// Loading data and normalization
	if( (spos+4)<nTimesamples ){
		ftemp1.x=d_input[gpos];	
		ftemp1.y=d_input[gpos+1];
		ftemp1.x = __fdividef( (ftemp1.x-signal_mean),signal_sd);	ftemp1.y = __fdividef( (ftemp1.y-signal_mean),signal_sd);
		ftemp2.x=d_input[gpos+2];
		ftemp2.y=d_input[gpos+3];
		ftemp2.x = __fdividef( (ftemp2.x-signal_mean),signal_sd);	ftemp2.y = __fdividef( (ftemp2.y-signal_mean),signal_sd);
		ftemp3.x=d_input[gpos+4];
		ftemp3.x = __fdividef( (ftemp3.x-signal_mean),signal_sd);

		Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
		Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
		Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
		Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
	}
	
	spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) &&  spos<dtm) d_decimated[blockIdx.y*dtm + spos]=Bw[1];
		
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];

	__syncthreads();
	
	taps1.x=1; taps1.y=1;
	ftemp1.x = Bw[0];
	ftemp1.y = Bw[4];
	SNR.x = ftemp1.x; SNR.y=ftemp1.y;
	taps.x = taps1.x; taps.y = taps1.y;
	
	Taps = 2;
	ftemp2 = Calculate_SNR_const_memory(Bw[1], Bw[5], Taps);
	Compare( &ftemp1, &taps1, &ftemp2, &Taps);
	
	// Using constant memory
	taps3.x = 3; taps3.y = 3; Taps = 4;
	ftemp3 = Calculate_SNR_const_memory(Bw[2], Bw[6], taps3.x);
	ftemp4 = Calculate_SNR_const_memory(Bw[3], Bw[7], Taps);
	Compare( &ftemp3, &taps3, &ftemp4, &Taps);
	
	Compare( &ftemp1, &taps1, &ftemp3, &taps3);
	Compare( &SNR, &taps, &ftemp1, &taps1);

	s_SNRs[threadIdx.x]=SNR;
	s_taps[threadIdx.x]=taps;
	
	for(d=4; d<nBoxcars; d=d+4){
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
			
			// Using constant memory
			taps1.x = d+1; taps1.y = d+1; Taps = d+2;
			ftemp1 = Calculate_SNR_const_memory(Bw[0], Bw[4], taps1.x);
			ftemp2 = Calculate_SNR_const_memory(Bw[1], Bw[5], Taps);
			Compare( &ftemp1, &taps1, &ftemp2, &Taps);

			// Using constant memory
			taps3.x = d+3; taps3.y = d+3; Taps = d+4;
			ftemp3 = Calculate_SNR_const_memory(Bw[2], Bw[6], taps3.x);
			ftemp4 = Calculate_SNR_const_memory(Bw[3], Bw[7], Taps);
			Compare( &ftemp3, &taps3, &ftemp4, &Taps);
			
			Compare( &ftemp1, &taps1, &ftemp3, &taps3);
			Compare( &SNR, &taps, &ftemp1, &taps1);
			
			s_SNRs[spos]=SNR;
			s_taps[spos]=taps;
		}
	}
	
	if(spos>=0) s_input[spos].x = Bw[3];
	
	__syncthreads();
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) && spos<nTimesamples){
		gpos=blockIdx.y*nTimesamples + spos;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = s_taps[threadIdx.x].x;
		d_output_taps[gpos + 1] = s_taps[threadIdx.x].y;
		
		spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
		if(spos<dtm) d_bv_out[blockIdx.y*dtm + spos] = s_input[threadIdx.x].x;
	}
}


__global__ void PD_GPU_Nth_BLN(float const* __restrict__ d_input, float *d_bv_in, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, const int nTimesamples, const int nBoxcars, const int startTaps, const int DIT_value, const int dtm) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_BV[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, ftemp4, SNR;
	ushort2 taps, taps1, taps3;
	ushort Taps;
	
	Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0;
	s_BV[threadIdx.x].x=0; s_BV[threadIdx.x].y=0;
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	gpos=blockIdx.y*nTimesamples + spos;
	
	if( (spos+4)<nTimesamples ){
		ftemp1.x=d_input[gpos];
		ftemp1.y=d_input[gpos + 1];
		ftemp2.x=d_input[gpos + 2];
		ftemp2.y=d_input[gpos + 3];
		ftemp3.x=d_input[gpos + 4];
		
		Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
		Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
		Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
		Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
		
		s_BV[threadIdx.x].x=d_bv_in[gpos];	
		s_BV[threadIdx.x].y=d_bv_in[gpos+1];
	}
	
	spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) &&  spos<dtm) d_decimated[blockIdx.y*dtm + spos]=Bw[1];
	
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];
	
	SNR.x=0; SNR.y=0;
	__syncthreads();
	
	taps1.x = startTaps + DIT_value; taps1.y = startTaps + DIT_value; Taps = startTaps + 2*DIT_value;
	ftemp1 = Calculate_SNR(s_BV[threadIdx.x].x + Bw[0], s_BV[threadIdx.x].y + Bw[4], taps1.x);
	SNR.x = ftemp1.x; SNR.y=ftemp1.y;
	taps.x = taps1.x; taps.y = taps1.y;
	ftemp2 = Calculate_SNR(s_BV[threadIdx.x].x + Bw[1], s_BV[threadIdx.x].y + Bw[5], Taps);
	Compare( &ftemp1, &taps1, &ftemp2, &Taps);

	taps3.x=startTaps + 3*DIT_value; taps3.y=startTaps + 3*DIT_value; Taps = startTaps + 4*DIT_value;
	ftemp3 = Calculate_SNR(s_BV[threadIdx.x].x + Bw[2], s_BV[threadIdx.x].y + Bw[6], taps3.x);
	ftemp4 = Calculate_SNR(s_BV[threadIdx.x].x + Bw[3], s_BV[threadIdx.x].y + Bw[7], Taps);
	Compare( &ftemp3, &taps3, &ftemp4, &Taps);
	
	Compare( &ftemp1, &taps1, &ftemp3, &taps3);
	Compare( &SNR, &taps, &ftemp1, &taps1);
	
	s_SNRs[threadIdx.x]=SNR;
	s_taps[threadIdx.x]=taps;
	
	for(d=4; d<nBoxcars; d=d+4){
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
			
			taps1.x = startTaps + d*DIT_value + DIT_value; taps1.y = startTaps + d*DIT_value + DIT_value; 
			Taps = startTaps + d*DIT_value + 2*DIT_value;
			ftemp1 = Calculate_SNR(s_BV[spos].x + Bw[0], s_BV[spos].y + Bw[4], taps1.x);
			ftemp2 = Calculate_SNR(s_BV[spos].x + Bw[1], s_BV[spos].y + Bw[5], Taps);
			Compare( &ftemp1, &taps1, &ftemp2, &Taps);
			
			taps3.x = startTaps + d*DIT_value + 3*DIT_value; taps3.y = startTaps + d*DIT_value + 3*DIT_value; 
			Taps = startTaps + d*DIT_value + 4*DIT_value;
			ftemp3 = Calculate_SNR(s_BV[spos].x + Bw[2], s_BV[spos].y + Bw[6], taps3.x);
			ftemp4 = Calculate_SNR(s_BV[spos].x + Bw[3], s_BV[spos].y + Bw[7], Taps);
			Compare( &ftemp3, &taps3, &ftemp4, &Taps);
			
			Compare( &ftemp1, &taps1, &ftemp3, &taps3);
			Compare( &SNR, &taps, &ftemp1, &taps1);
			
			s_SNRs[spos]=SNR;
			s_taps[spos]=taps;
		}

	}
	
	if(spos>=0) s_input[spos].x = s_BV[spos].x + Bw[3];
	
	__syncthreads();
	
	spos = blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1))  && spos<(nTimesamples) ){
		gpos=blockIdx.y*nTimesamples + spos;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = s_taps[threadIdx.x].x;
		d_output_taps[gpos + 1] = s_taps[threadIdx.x].y;
		
		spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
		if(spos<dtm) d_bv_out[blockIdx.y*dtm + spos] = s_input[threadIdx.x].x;
	}
}


__global__ void PD_GPU_Nth_BLN_EACH(float const* __restrict__ d_input, float *d_bv_in, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD_BV, float *d_MSD_DIT, const int nTimesamples, const int nBoxcars, const int startTaps, const int DIT_value, const int dtm) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_BV[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	float mean_BV  = d_MSD_BV[0];
	float mean_DIT = d_MSD_DIT[0];
	float sd_BV    = d_MSD_BV[1]*d_MSD_BV[1];
	float sd_DIT   = d_MSD_DIT[1]*d_MSD_DIT[1];
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, ftemp4, SNR;
	ushort2 taps, taps1, taps3;
	ushort Taps;
	
	Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0;
	s_BV[threadIdx.x].x=0; s_BV[threadIdx.x].y=0;
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	gpos=blockIdx.y*nTimesamples + spos;
	
	if( (spos+4)<nTimesamples ){
		ftemp1.x=d_input[gpos];	
		ftemp1.y=d_input[gpos+1];
		ftemp2.x=d_input[gpos+2];
		ftemp2.y=d_input[gpos+3];
		ftemp3.x=d_input[gpos+4];
				
		
		Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
		Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
		Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
		Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
		
		s_BV[threadIdx.x].x=d_bv_in[gpos];	
		s_BV[threadIdx.x].y=d_bv_in[gpos+1];
	}
	
	spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) &&  spos<dtm) d_decimated[blockIdx.y*dtm + spos]=Bw[1];
	
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];
	
	SNR.x=0; SNR.y=0;
	__syncthreads();
	
	//----------- Old
	//taps1.x = 1; taps1.y=2;
	//ftemp1.x = __frsqrt_rn(taps1.x*sd_DIT)*(s_BV[threadIdx.x].x + Bw[0]);     ftemp1.y = __frsqrt_rn(taps1.x*sd_DIT)*(s_BV[threadIdx.x].y + Bw[4]);
	//ftemp2.x = __frsqrt_rn(taps1.y*sd_DIT)*(s_BV[threadIdx.x].x + Bw[1]);     ftemp2.y = __frsqrt_rn(taps1.y*sd_DIT)*(s_BV[threadIdx.x].y + Bw[5]);
	//if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=2;} else taps1.x=1;
	//if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=2;} else taps1.y=1;
	//
	//taps2.x=3; taps2.y=4;
	//ftemp2.x = __frsqrt_rn(taps2.x*sd_DIT)*(s_BV[threadIdx.x].x + Bw[2]);     ftemp2.y = __frsqrt_rn(taps2.x*sd_DIT)*(s_BV[threadIdx.x].y + Bw[6]);
	//ftemp3.x = __frsqrt_rn(taps2.y*sd_DIT)*(s_BV[threadIdx.x].x + Bw[3]);     ftemp3.y = __frsqrt_rn(taps2.y*sd_DIT)*(s_BV[threadIdx.x].y + Bw[7]);
	//if(ftemp3.x>ftemp2.x) {ftemp2.x=ftemp3.x; taps2.x=4;} else taps2.x=3;
	//if(ftemp3.y>ftemp2.y) {ftemp2.y=ftemp3.y; taps2.y=4;} else taps2.y=3;
	//---------- Old
	
	taps1.x = 1; taps1.y = 1; Taps = 2;
	ftemp1 = Calculate_SNR_EACH( (s_BV[threadIdx.x].x + Bw[0])-mean_BV-taps1.x*mean_DIT, (s_BV[threadIdx.x].y + Bw[4])-mean_BV-taps1.x*mean_DIT, sd_BV + taps1.x*sd_DIT);
	SNR.x = ftemp1.x; SNR.y=ftemp1.y;
	taps.x = taps1.x; taps.y = taps1.y;
	ftemp2 = Calculate_SNR_EACH( (s_BV[threadIdx.x].x + Bw[1])-mean_BV-Taps*mean_DIT, (s_BV[threadIdx.x].y + Bw[5])-mean_BV-Taps*mean_DIT, sd_BV + Taps*sd_DIT);
	Compare( &ftemp1, &taps1, &ftemp2, &Taps);

	taps3.x=3; taps3.y=3; Taps = 4;
	ftemp3 = Calculate_SNR_EACH( (s_BV[threadIdx.x].x + Bw[2])-mean_BV-taps3.x*mean_DIT, (s_BV[threadIdx.x].y + Bw[6])-mean_BV-taps3.x*mean_DIT, sd_BV + taps3.x*sd_DIT);
	ftemp4 = Calculate_SNR_EACH( (s_BV[threadIdx.x].x + Bw[3])-mean_BV-Taps*mean_DIT, (s_BV[threadIdx.x].y + Bw[7])-mean_BV-Taps*mean_DIT, sd_BV + Taps*sd_DIT);
	Compare( &ftemp3, &taps3, &ftemp4, &Taps);
	
	Compare( &ftemp1, &taps1, &ftemp3, &taps3);
	Compare( &SNR, &taps, &ftemp1, &taps1);
	
	s_SNRs[threadIdx.x]=SNR;
	s_taps[threadIdx.x]=taps;
	
	for(d=4; d<nBoxcars; d=d+4){
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
			
			//------------- old
			//Taps = d;
			//taps1.x = Taps + 1; taps1.y=Taps + 2;
			//ftemp1.x = __frsqrt_rn(taps1.x*sd_DIT)*(s_BV[spos].x + Bw[0]);     ftemp1.y = __frsqrt_rn(taps1.x*sd_DIT)*(s_BV[spos].y + Bw[4]);
			//ftemp2.x = __frsqrt_rn(taps1.y*sd_DIT)*(s_BV[spos].x + Bw[1]);     ftemp2.y = __frsqrt_rn(taps1.y*sd_DIT)*(s_BV[spos].y + Bw[5]);
			//if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = Taps +2;} else taps1.x = Taps + 1;
			//if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = Taps +2;} else taps1.y = Taps + 1;
			//
			//taps2.x=Taps + 3; taps2.y=Taps + 4;
			//ftemp2.x = __frsqrt_rn(taps2.x*sd_DIT)*(s_BV[spos].x + Bw[2]);     ftemp2.y = __frsqrt_rn(taps2.x*sd_DIT)*(s_BV[spos].y + Bw[6]);
			//ftemp3.x = __frsqrt_rn(taps2.y*sd_DIT)*(s_BV[spos].x + Bw[3]);     ftemp3.y = __frsqrt_rn(taps2.y*sd_DIT)*(s_BV[spos].y + Bw[7]);
			//if(ftemp3.x>ftemp2.x) {ftemp2.x = ftemp3.x; taps2.x = Taps + 4;} else taps2.x = Taps + 3;
			//if(ftemp3.y>ftemp2.y) {ftemp2.y = ftemp3.y; taps2.y = Taps + 4;} else taps2.y = Taps + 3;
			//
			//if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = taps2.x;}
			//if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = taps2.y;}
			//
			//if(ftemp1.x>SNR.x) {SNR.x = ftemp1.x; taps.x = taps1.x;}
			//if(ftemp1.y>SNR.y) {SNR.y = ftemp1.y; taps.y = taps1.y;}
			//--------------- old
			
			taps1.x = d+1; taps1.y = d+1; 
			Taps = d+2;
			ftemp1 = Calculate_SNR_EACH( (s_BV[spos].x + Bw[0])-mean_BV-taps1.x*mean_DIT, (s_BV[spos].y + Bw[4])-mean_BV-taps1.x*mean_DIT, sd_BV + taps1.x*sd_DIT);
			ftemp2 = Calculate_SNR_EACH( (s_BV[spos].x + Bw[1])-mean_BV-Taps*mean_DIT, (s_BV[spos].y + Bw[5])-mean_BV-Taps*mean_DIT, sd_BV + Taps*sd_DIT);
			Compare( &ftemp1, &taps1, &ftemp2, &Taps);
			
			taps3.x = d+3; taps3.y = d+3; 
			Taps = d+4;
			ftemp3 = Calculate_SNR_EACH( (s_BV[spos].x + Bw[2])-mean_BV-taps3.x*mean_DIT, (s_BV[spos].y + Bw[6])-mean_BV-taps3.x*mean_DIT, sd_BV + taps3.x*sd_DIT);
			ftemp4 = Calculate_SNR_EACH( (s_BV[spos].x + Bw[3])-mean_BV-Taps*mean_DIT, (s_BV[spos].y + Bw[7])-mean_BV-Taps*mean_DIT, sd_BV + Taps*sd_DIT);
			Compare( &ftemp3, &taps3, &ftemp4, &Taps);
			
			Compare( &ftemp1, &taps1, &ftemp3, &taps3);
			Compare( &SNR, &taps, &ftemp1, &taps1);
			
			s_SNRs[spos]=SNR;
			s_taps[spos]=taps;
		}

	}
	
	if(spos>=0) s_input[spos].x = s_BV[spos].x + Bw[3];
	
	__syncthreads();
	
	spos = blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1))  && spos<(nTimesamples) ){
		gpos=blockIdx.y*nTimesamples + spos;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = startTaps + s_taps[threadIdx.x].x*DIT_value;
		d_output_taps[gpos + 1] = startTaps + s_taps[threadIdx.x].y*DIT_value;
		
		spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
		if(spos<dtm) d_bv_out[blockIdx.y*dtm + spos] = s_input[threadIdx.x].x;
	}
}

//*********************************** BLN *********************************************
//*************************************************************************************




// *************************************************************************************************************
// *************************************************************************************************************
// *************************************************************************************************************

__global__ void PD_GPU_1st_LA(float const* __restrict__ d_input, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int nTimesamples, int nBoxcars, const int dtm) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, SNR;
	ushort2 taps, taps1, taps2;
	float signal_mean=d_MSD[0];
	float signal_sd=d_MSD[1];
	float modifier=d_MSD[2];
	
	Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0; SNR.x=0; SNR.y=0;
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	gpos = blockIdx.y*nTimesamples + spos;
	
	// Loading data and normalization
	if( (spos+4)<nTimesamples ){
		ftemp1.x=d_input[gpos];	
		ftemp1.y=d_input[gpos+1];
		ftemp2.x=d_input[gpos+2];
		ftemp2.y=d_input[gpos+3];
		ftemp3.x=d_input[gpos+4];

		Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
		Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
		Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
		Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
	}
	
	spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) &&  spos<dtm) d_decimated[blockIdx.y*dtm + spos]=Bw[1];
		
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];

	__syncthreads();
	
	ftemp1.x = __fdividef(Bw[0]-signal_mean, signal_sd);                 ftemp1.y = __fdividef(Bw[4]-signal_mean, signal_sd);
	ftemp2.x = __fdividef(Bw[1]-2*signal_mean,signal_sd + modifier);     ftemp2.y = __fdividef(Bw[5]-2*signal_mean,signal_sd + modifier);
	
	if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=2;} else taps1.x=1;
	if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=2;} else taps1.y=1;
	

	ftemp2.x = __fdividef(Bw[2]-3*signal_mean, signal_sd + 2*modifier);    ftemp2.y = __fdividef(Bw[6]-3*signal_mean, signal_sd + 2*modifier);
	ftemp3.x = __fdividef(Bw[3]-4*signal_mean, signal_sd + 3*modifier);    ftemp3.y = __fdividef(Bw[7]-4*signal_mean, signal_sd + 3*modifier);
	
	if(ftemp3.x>ftemp2.x) {ftemp2.x=ftemp3.x; taps2.x=4;} else taps2.x=3;
	if(ftemp3.y>ftemp2.y) {ftemp2.y=ftemp3.y; taps2.y=4;} else taps2.y=3;
	
	if(ftemp1.x>ftemp2.x) {SNR.x=ftemp1.x; taps.x=taps1.x;} else {SNR.x=ftemp2.x; taps.x=taps2.x;}
	if(ftemp1.y>ftemp2.y) {SNR.y=ftemp1.y; taps.y=taps1.y;} else {SNR.y=ftemp2.y; taps.y=taps2.y;}

	s_SNRs[threadIdx.x]=SNR;
	s_taps[threadIdx.x]=taps;
	
	for(d=4; d<nBoxcars; d=d+4){
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
			
			// Using constant memory
			ftemp1.x = __fdividef( (Bw[0]-(d+1)*signal_mean) ,(signal_sd + d*modifier));
			ftemp1.y = __fdividef( (Bw[4]-(d+1)*signal_mean) ,(signal_sd + d*modifier));
			ftemp2.x = __fdividef( (Bw[1]-(d+2)*signal_mean) ,(signal_sd + (d+1)*modifier));
			ftemp2.y = __fdividef( (Bw[5]-(d+2)*signal_mean) ,(signal_sd + (d+1)*modifier));
			if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = d+2;} else taps1.x = d+1;
			if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = d+2;} else taps1.y = d+1;

			// Using constant memory
			ftemp2.x = __fdividef( (Bw[2]-(d+3)*signal_mean) ,(signal_sd + (d+2)*modifier));
			ftemp2.y = __fdividef( (Bw[6]-(d+3)*signal_mean),(signal_sd + (d+2)*modifier));
			ftemp3.x = __fdividef( (Bw[3]-(d+4)*signal_mean) ,(signal_sd + (d+3)*modifier));
			ftemp3.y = __fdividef( (Bw[7]-(d+4)*signal_mean),(signal_sd + (d+3)*modifier));
			if(ftemp3.x>ftemp2.x) {ftemp2.x = ftemp3.x; taps2.x = d+4;} else taps2.x = d+3;
			if(ftemp3.y>ftemp2.y) {ftemp2.y = ftemp3.y; taps2.y = d+4;} else taps2.y = d+3;
			
			
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
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) && spos<nTimesamples){
		gpos=blockIdx.y*nTimesamples + spos;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = s_taps[threadIdx.x].x;
		d_output_taps[gpos + 1] = s_taps[threadIdx.x].y;
		
		spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
		if(spos<dtm) d_bv_out[blockIdx.y*dtm + spos] = s_input[threadIdx.x].x;
	}
}



__global__ void PD_GPU_Nth_LA(float const* __restrict__ d_input, float *d_bv_in, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, const int nTimesamples, const int nBoxcars, const int startTaps, const int DIT_value, const int dtm) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_BV[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, SNR;
	ushort2 taps, taps1, taps2;
	int Taps;
	float signal_mean=d_MSD[0];
	float signal_sd=d_MSD[1];
	float modifier=d_MSD[2];
	
	Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0;
	s_BV[threadIdx.x].x=0; s_BV[threadIdx.x].y=0;
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	gpos=blockIdx.y*nTimesamples + spos;
	
	if( (spos+4)<nTimesamples ){
		ftemp1.x=d_input[gpos];
		ftemp1.y=d_input[gpos + 1];
		ftemp2.x=d_input[gpos + 2];
		ftemp2.y=d_input[gpos + 3];
		ftemp3.x=d_input[gpos + 4];
		
		
		Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
		Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
		Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
		Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
		
		s_BV[threadIdx.x].x=d_bv_in[gpos];	
		s_BV[threadIdx.x].y=d_bv_in[gpos+1];
	}
	
	spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) &&  spos<dtm) d_decimated[blockIdx.y*dtm + spos]=Bw[1];

	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];
	
	SNR.x=0; SNR.y=0;
	__syncthreads();
	
	taps1.x = startTaps + DIT_value; taps1.y=startTaps + 2*DIT_value;
	ftemp1.x = __fdividef( ((s_BV[threadIdx.x].x + Bw[0]) - taps1.x*signal_mean) , (signal_sd + (taps1.x - 1)*modifier));       ftemp1.y = __fdividef( ((s_BV[threadIdx.x].y + Bw[4]) - taps1.x*signal_mean) , (signal_sd + (taps1.x - 1)*modifier));
	ftemp2.x = __fdividef( ((s_BV[threadIdx.x].x + Bw[1]) - taps1.y*signal_mean) , (signal_sd + (taps1.y - 1)*modifier));       ftemp2.y = __fdividef( ((s_BV[threadIdx.x].y + Bw[5]) - taps1.y*signal_mean) , (signal_sd + (taps1.y - 1)*modifier));
	if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=startTaps + 2*DIT_value;} else taps1.x=startTaps + DIT_value;
	if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=startTaps + 2*DIT_value;} else taps1.y=startTaps + DIT_value;

	taps2.x=startTaps + 3*DIT_value; taps2.y=startTaps + 4*DIT_value;
	ftemp2.x = __fdividef( ((s_BV[threadIdx.x].x + Bw[2]) - taps2.x*signal_mean) , (signal_sd + (taps2.x - 1)*modifier));       ftemp2.y = __fdividef( ((s_BV[threadIdx.x].y + Bw[6]) - taps2.x*signal_mean) , (signal_sd + (taps2.x - 1)*modifier));
	ftemp3.x = __fdividef( ((s_BV[threadIdx.x].x + Bw[3]) - taps2.y*signal_mean) , (signal_sd + (taps2.y - 1)*modifier));       ftemp3.y = __fdividef( ((s_BV[threadIdx.x].y + Bw[7]) - taps2.y*signal_mean) , (signal_sd + (taps2.y - 1)*modifier));
	if(ftemp3.x>ftemp2.x) {ftemp2.x=ftemp3.x; taps2.x=startTaps+4*DIT_value;} else taps2.x=startTaps+3*DIT_value;
	if(ftemp3.y>ftemp2.y) {ftemp2.y=ftemp3.y; taps2.y=startTaps+4*DIT_value;} else taps2.y=startTaps+3*DIT_value;
	
	if(ftemp1.x>ftemp2.x) {SNR.x=ftemp1.x; taps.x=taps1.x;} else {SNR.x=ftemp2.x; taps.x=taps2.x;}
	if(ftemp1.y>ftemp2.y) {SNR.y=ftemp1.y; taps.y=taps1.y;} else {SNR.y=ftemp2.y; taps.y=taps2.y;}
	
	s_SNRs[threadIdx.x]=SNR;
	s_taps[threadIdx.x]=taps;
	
	for(d=4; d<nBoxcars; d=d+4){
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
			
			Taps=startTaps + d*DIT_value;
			
			taps1.x = Taps + DIT_value; taps1.y=Taps + 2*DIT_value;
			ftemp1.x = __fdividef( ((s_BV[spos].x + Bw[0]) - taps1.x*signal_mean) ,(signal_sd + (taps1.x - 1)*modifier));      ftemp1.y = __fdividef( ((s_BV[spos].y + Bw[4]) - taps1.x*signal_mean) ,(signal_sd + (taps1.x - 1)*modifier));
			ftemp2.x = __fdividef( ((s_BV[spos].x + Bw[1]) - taps1.y*signal_mean) ,(signal_sd + (taps1.y - 1)*modifier));      ftemp2.y = __fdividef( ((s_BV[spos].y + Bw[5]) - taps1.y*signal_mean) ,(signal_sd + (taps1.y - 1)*modifier));
			if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = Taps +2*DIT_value;} else taps1.x = Taps + DIT_value;
			if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = Taps +2*DIT_value;} else taps1.y = Taps + DIT_value;
			
			taps2.x=Taps + 3*DIT_value; taps2.y=Taps + 4*DIT_value;
			ftemp2.x = __fdividef( ((s_BV[spos].x + Bw[2]) - taps2.x*signal_mean) ,(signal_sd + (taps2.x - 1)*modifier));      ftemp2.y = __fdividef( ((s_BV[spos].y + Bw[6]) - taps2.x*signal_mean) ,(signal_sd + (taps2.x - 1)*modifier));
			ftemp3.x = __fdividef( ((s_BV[spos].x + Bw[3]) - taps2.y*signal_mean) ,(signal_sd + (taps2.y - 1)*modifier));      ftemp3.y = __fdividef( ((s_BV[spos].y + Bw[7]) - taps2.y*signal_mean) ,(signal_sd + (taps2.y - 1)*modifier));
			if(ftemp3.x>ftemp2.x) {ftemp2.x = ftemp3.x; taps2.x = Taps + 4*DIT_value;} else taps2.x = Taps + 3*DIT_value;
			if(ftemp3.y>ftemp2.y) {ftemp2.y = ftemp3.y; taps2.y = Taps + 4*DIT_value;} else taps2.y = Taps + 3*DIT_value;
			
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
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1))  && spos<(nTimesamples) ){
		gpos=blockIdx.y*nTimesamples + spos;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = s_taps[threadIdx.x].x;
		d_output_taps[gpos + 1] = s_taps[threadIdx.x].y;
		
		spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
		if(spos<dtm) d_bv_out[blockIdx.y*dtm + spos] = s_input[threadIdx.x].x;
	}
}



__global__ void PD_GPU_Nth_LA_EACH(float const* __restrict__ d_input, float *d_bv_in, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, const int nTimesamples, const int nBoxcars, const int startTaps, const int DIT_value, const int dtm) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_BV[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, SNR;
	ushort2 taps, taps1, taps2;
	float BV_mean=d_MSD[0];
	float DIT_mean=d_MSD[3];
	float signal_sd=d_MSD[1];
	float modifier=d_MSD[2];
	
	Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0;
	s_BV[threadIdx.x].x=0; s_BV[threadIdx.x].y=0;
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	gpos=blockIdx.y*nTimesamples + spos;
	
	if( (spos+4)<nTimesamples ){
		ftemp1.x=d_input[gpos];
		ftemp1.y=d_input[gpos + 1];
		ftemp2.x=d_input[gpos + 2];
		ftemp2.y=d_input[gpos + 3];
		ftemp3.x=d_input[gpos + 4];
		
		
		Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
		Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
		Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
		Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
		
		s_BV[threadIdx.x].x=d_bv_in[gpos]   - BV_mean;	
		s_BV[threadIdx.x].y=d_bv_in[gpos+1] - BV_mean;
	}
	
	spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) &&  spos<dtm) d_decimated[blockIdx.y*dtm + spos]=Bw[1];
	
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];
	
	SNR.x=0; SNR.y=0;
	__syncthreads();
	
	taps1.x = 1; taps1.y=2;
	ftemp1.x = __fdividef( (s_BV[threadIdx.x].x + (Bw[0] - taps1.x*DIT_mean)) , (signal_sd + taps1.x*modifier));
	ftemp1.y = __fdividef( (s_BV[threadIdx.x].y + (Bw[4] - taps1.x*DIT_mean)) , (signal_sd + taps1.x*modifier));
	ftemp2.x = __fdividef( (s_BV[threadIdx.x].x + (Bw[1] - taps1.y*DIT_mean)) , (signal_sd + taps1.y*modifier));
	ftemp2.y = __fdividef( (s_BV[threadIdx.x].y + (Bw[5] - taps1.y*DIT_mean)) , (signal_sd + taps1.y*modifier));
	if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=2;} else taps1.x=1;
	if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=2;} else taps1.y=1;

	taps2.x=3; taps2.y=4;
	ftemp2.x = __fdividef( (s_BV[threadIdx.x].x + (Bw[2] - taps2.x*DIT_mean)) , (signal_sd + taps2.x*modifier));
	ftemp2.y = __fdividef( (s_BV[threadIdx.x].y + (Bw[6] - taps2.x*DIT_mean)) , (signal_sd + taps2.x*modifier));
	ftemp3.x = __fdividef( (s_BV[threadIdx.x].x + (Bw[3] - taps2.y*DIT_mean)) , (signal_sd + taps2.y*modifier));
	ftemp3.y = __fdividef( (s_BV[threadIdx.x].y + (Bw[7] - taps2.y*DIT_mean)) , (signal_sd + taps2.y*modifier));
	if(ftemp3.x>ftemp2.x) {ftemp2.x=ftemp3.x; taps2.x=4;} else taps2.x=3;
	if(ftemp3.y>ftemp2.y) {ftemp2.y=ftemp3.y; taps2.y=4;} else taps2.y=3;
	
	if(ftemp1.x>ftemp2.x) {SNR.x=ftemp1.x; taps.x=taps1.x;} else {SNR.x=ftemp2.x; taps.x=taps2.x;}
	if(ftemp1.y>ftemp2.y) {SNR.y=ftemp1.y; taps.y=taps1.y;} else {SNR.y=ftemp2.y; taps.y=taps2.y;}
	
	s_SNRs[threadIdx.x]=SNR;
	s_taps[threadIdx.x]=taps;
	
	for(d=4; d<nBoxcars; d=d+4){
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
			
			taps1.x = d+1; taps1.y=d+2;
			ftemp1.x = __fdividef( (s_BV[spos].x + (Bw[0] - taps1.x*DIT_mean)) ,(signal_sd + taps1.x*modifier));      ftemp1.y = __fdividef( (s_BV[spos].y + (Bw[4] - taps1.x*DIT_mean)) ,(signal_sd + taps1.x*modifier));
			ftemp2.x = __fdividef( (s_BV[spos].x + (Bw[1] - taps1.y*DIT_mean)) ,(signal_sd + taps1.y*modifier));      ftemp2.y = __fdividef( (s_BV[spos].y + (Bw[5] - taps1.y*DIT_mean)) ,(signal_sd + taps1.y*modifier));
			if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = d+2;} else taps1.x = d+1;
			if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = d+2;} else taps1.y = d+1;
			
			taps2.x=d+3; taps2.y=d+4;
			ftemp2.x = __fdividef( (s_BV[spos].x + (Bw[2] - taps2.x*DIT_mean)) ,(signal_sd + taps2.x*modifier));      ftemp2.y = __fdividef( (s_BV[spos].y + (Bw[6] - taps2.x*DIT_mean)) ,(signal_sd + taps2.x*modifier));
			ftemp3.x = __fdividef( (s_BV[spos].x + (Bw[3] - taps2.y*DIT_mean)) ,(signal_sd + taps2.y*modifier));      ftemp3.y = __fdividef( (s_BV[spos].y + (Bw[7] - taps2.y*DIT_mean)) ,(signal_sd + taps2.y*modifier));
			if(ftemp3.x>ftemp2.x) {ftemp2.x = ftemp3.x; taps2.x = d+4;} else taps2.x = d+3;
			if(ftemp3.y>ftemp2.y) {ftemp2.y = ftemp3.y; taps2.y = d+4;} else taps2.y = d+3;
			
			if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = taps2.x;}
			if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = taps2.y;}
			
			if(ftemp1.x>SNR.x) {SNR.x = ftemp1.x; taps.x = taps1.x;}
			if(ftemp1.y>SNR.y) {SNR.y = ftemp1.y; taps.y = taps1.y;}
			
			s_SNRs[spos]=SNR;
			s_taps[spos]=taps;
		}

	}
	
	if(spos>=0) s_input[spos].x = s_BV[spos].x + BV_mean + Bw[3];
	
	__syncthreads();
	
	spos = blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1))  && spos<(nTimesamples) ){
		gpos=blockIdx.y*nTimesamples + spos;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = startTaps + s_taps[threadIdx.x].x*DIT_value;
		d_output_taps[gpos + 1] = startTaps + s_taps[threadIdx.x].y*DIT_value;
		
		spos = blockIdx.x*(PD_NTHREADS-(nBoxcars>>1)) + threadIdx.x;
		if(spos<dtm) d_bv_out[blockIdx.y*dtm + spos] = s_input[threadIdx.x].x;
	}
}



#endif
