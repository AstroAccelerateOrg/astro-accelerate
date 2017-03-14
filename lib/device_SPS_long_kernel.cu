// Added by Karel Adamek 

#ifndef SPS_LONG_KERNEL_H_
#define SPS_LONG_KERNEL_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include "headers/params.h"

__global__ void PD_GPU_1st_float1(float const* __restrict__ d_input, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int nTimesamples, int nBoxcars) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ int2 s_taps[PD_NTHREADS];
	
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, SNR;
	int2 taps, taps1, taps2;
	float signal_mean=d_MSD[0];
	float signal_sd=d_MSD[2];
	float modifier=d_MSD[1];
	
	gpos=blockIdx.y*nTimesamples + blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	
	ftemp1.x=d_input[gpos];
	ftemp1.y=d_input[gpos+1];
	ftemp2.x=d_input[gpos+2];
	ftemp2.y=d_input[gpos+3];
	ftemp3.x=d_input[gpos+4];
	
	Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
	Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
	Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
	Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
	
	d_decimated[(gpos>>1)]=Bw[1];
	
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];

	SNR.x=0; SNR.y=0;
	__syncthreads();
	
	ftemp1.x = __fdividef( (Bw[0]-signal_mean),signal_sd);                 ftemp1.y = __fdividef( (Bw[4]-signal_mean),signal_sd);
	ftemp2.x = __fdividef( (Bw[1]-2*signal_mean),(signal_sd + modifier));  ftemp2.y = __fdividef( (Bw[5]-2*signal_mean),(signal_sd + modifier));
	if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=2;} else taps1.x=1;
	if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=2;} else taps1.y=1;
	
	
	ftemp2.x = __fdividef( (Bw[2]-3*signal_mean),(signal_sd + 2*modifier)); ftemp2.y = __fdividef( (Bw[6]-3*signal_mean),(signal_sd + 2*modifier));
	ftemp3.x = __fdividef( (Bw[3]-4*signal_mean),(signal_sd + 3*modifier)); ftemp3.y = __fdividef( (Bw[7]-4*signal_mean),(signal_sd + 3*modifier));
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
			
			
			ftemp1.x = __fdividef( (Bw[0]-(d+1)*signal_mean) ,(signal_sd + d*modifier));     ftemp1.y = __fdividef( (Bw[4]-(d+1)*signal_mean) ,(signal_sd + d*modifier));
			ftemp2.x = __fdividef( (Bw[1]-(d+2)*signal_mean) ,(signal_sd + (d+1)*modifier)); ftemp2.y = __fdividef( (Bw[5]-(d+2)*signal_mean) ,(signal_sd + (d+1)*modifier));
			if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = d+2;} else taps1.x = d+1;
			if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = d+2;} else taps1.y = d+1;

			
			ftemp2.x = __fdividef( (Bw[2]-(d+3)*signal_mean) ,(signal_sd + (d+2)*modifier)); ftemp2.y = __fdividef( (Bw[6]-(d+3)*signal_mean),(signal_sd + (d+2)*modifier));
			ftemp3.x = __fdividef( (Bw[3]-(d+4)*signal_mean) ,(signal_sd + (d+3)*modifier)); ftemp3.y = __fdividef( (Bw[7]-(d+4)*signal_mean),(signal_sd + (d+3)*modifier));
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
	
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) ){
		gpos=blockIdx.y*nTimesamples + blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = (ushort) s_taps[threadIdx.x].x;
		d_output_taps[gpos + 1] = (ushort) s_taps[threadIdx.x].y;
		
		d_bv_out[(gpos>>1)] = s_input[threadIdx.x].x;
	}
}



__global__ void PD_GPU_Nth_float1(float const* __restrict__ d_input, float *d_bv_in, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int nTimesamples, int nBoxcars, int startTaps, int DIT_value) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_BV[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ int2 s_taps[PD_NTHREADS];
	
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, SNR;
	int2 taps, taps1, taps2;
	int Taps;
	float signal_mean=d_MSD[0];
	float signal_sd=d_MSD[2];
	float modifier=d_MSD[1];
	
	gpos=blockIdx.y*nTimesamples + blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	
	ftemp1.x=d_input[gpos];
	ftemp1.y=d_input[gpos + 1];
	ftemp2.x=d_input[gpos + 2];
	ftemp2.y=d_input[gpos + 3];
	ftemp3.x=d_input[gpos + 4];
	
	Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
	Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
	Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
	Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
	
	d_decimated[(gpos>>1)]=Bw[1];

	
	s_BV[threadIdx.x].x=d_bv_in[gpos];	
	s_BV[threadIdx.x].y=d_bv_in[gpos+1];
	
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
	ftemp2.x = __fdividef( ((s_BV[threadIdx.x].x + Bw[2]) - taps2.x*signal_mean) ,(signal_sd + (taps2.x - 1)*modifier)); ftemp2.y = __fdividef( ((s_BV[threadIdx.x].y + Bw[6]) - taps2.x*signal_mean) ,(signal_sd + (taps2.x - 1)*modifier));
	ftemp3.x = __fdividef( ((s_BV[threadIdx.x].x + Bw[3]) - taps2.y*signal_mean) ,(signal_sd + (taps2.y - 1)*modifier)); ftemp3.y = __fdividef( ((s_BV[threadIdx.x].y + Bw[7]) - taps2.y*signal_mean) ,(signal_sd + (taps2.y - 1)*modifier));
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
	
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) ){
		gpos=blockIdx.y*nTimesamples + blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = (ushort) s_taps[threadIdx.x].x;
		d_output_taps[gpos + 1] = (ushort) s_taps[threadIdx.x].y;
		
		d_bv_out[(gpos>>1)] = s_input[threadIdx.x].x;
	}
}



//*********************************************************************************************************************
//*********************************************************************************************************************
//*********************************************************************************************************************


__global__ void PD_GPU_1st_float1_BLN(float const* __restrict__ d_input, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int nTimesamples, int nBoxcars) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, SNR;
	ushort2 taps, taps1, taps2;
	float signal_mean=d_MSD[0];
	float signal_sd=d_MSD[1];
	
	gpos=blockIdx.y*nTimesamples + blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	
	// Loading data and normalization
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
	
	d_decimated[(gpos>>1)]=Bw[1];
	
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];

	SNR.x=0; SNR.y=0;
	__syncthreads();
	
	ftemp1.x = Bw[0];                    ftemp1.y = Bw[4];
	//ftemp2.x = __frsqrt_rn(2)*Bw[1];     ftemp2.y = __frsqrt_rn(2)*Bw[5];
	ftemp2.x = __fdividef(Bw[1],c_sqrt_taps[2]);     ftemp2.y = __fdividef(Bw[5],c_sqrt_taps[2]);
	
	if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=2;} else taps1.x=1;
	if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=2;} else taps1.y=1;
	
	// Using fast intrinsic
	//ftemp2.x = __frsqrt_rn(3)*Bw[2]; ftemp2.y = __frsqrt_rn(3)*Bw[6];
	//ftemp3.x = __frsqrt_rn(4)*Bw[3]; ftemp3.y = __frsqrt_rn(4)*Bw[7];
	// Using constant memory
	ftemp2.x = __fdividef(Bw[2],c_sqrt_taps[3]);    ftemp2.y = __fdividef(Bw[6],c_sqrt_taps[3]);
	ftemp3.x = __fdividef(Bw[3],c_sqrt_taps[4]);    ftemp3.y = __fdividef(Bw[7],c_sqrt_taps[4]);
	
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
			
			// Using fast intrinsic
			//ftemp1.x = __frsqrt_rn(d+1)*Bw[0];     ftemp1.y = __frsqrt_rn(d+1)*Bw[4];
			//ftemp2.x = __frsqrt_rn(d+2)*Bw[1];     ftemp2.y = __frsqrt_rn(d+2)*Bw[5];
			// Using constant memory
			ftemp1.x = __fdividef(Bw[0],c_sqrt_taps[d+1]);    ftemp1.y = __fdividef(Bw[4],c_sqrt_taps[d+1]);
			ftemp2.x = __fdividef(Bw[1],c_sqrt_taps[d+2]);    ftemp2.y = __fdividef(Bw[5],c_sqrt_taps[d+2]);
			
			if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = d+2;} else taps1.x = d+1;
			if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = d+2;} else taps1.y = d+1;

			// Using fast intrinsic
			//ftemp2.x = __frsqrt_rn(d+3)*Bw[2];     ftemp2.y = __frsqrt_rn(d+3)*Bw[6];
			//ftemp3.x = __frsqrt_rn(d+4)*Bw[3];     ftemp3.y = __frsqrt_rn(d+4)*Bw[7];
			// Using constant memory
			ftemp2.x = __fdividef(Bw[2],c_sqrt_taps[d+3]);    ftemp2.y = __fdividef(Bw[6],c_sqrt_taps[d+3]);
			ftemp3.x = __fdividef(Bw[3],c_sqrt_taps[d+4]);    ftemp3.y = __fdividef(Bw[7],c_sqrt_taps[d+4]);
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
	
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) ){
		gpos=blockIdx.y*nTimesamples + blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = s_taps[threadIdx.x].x;
		d_output_taps[gpos + 1] = s_taps[threadIdx.x].y;
		
		d_bv_out[(gpos>>1)] = s_input[threadIdx.x].x;
	}
}



__global__ void PD_GPU_Nth_float1_BLN(float const* __restrict__ d_input, float *d_bv_in, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, int nTimesamples, int nBoxcars, int startTaps, int DIT_value) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_BV[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, SNR;
	ushort2 taps, taps1, taps2;
	int Taps;
	
	gpos=blockIdx.y*nTimesamples + blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	
	ftemp1.x=d_input[gpos];
	ftemp1.y=d_input[gpos + 1];
	ftemp2.x=d_input[gpos + 2];
	ftemp2.y=d_input[gpos + 3];
	ftemp3.x=d_input[gpos + 4];
	
	Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
	Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
	Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
	Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
	
	d_decimated[(gpos>>1)]=Bw[1];

	
	s_BV[threadIdx.x].x=d_bv_in[gpos];	
	s_BV[threadIdx.x].y=d_bv_in[gpos+1];
	
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];
	
	SNR.x=0; SNR.y=0;
	__syncthreads();
	
	taps1.x = startTaps + DIT_value; taps1.y=startTaps + 2*DIT_value;
	ftemp1.x = __frsqrt_rn(taps1.x)*(s_BV[threadIdx.x].x + Bw[0]);     ftemp1.y = __frsqrt_rn(taps1.x)*(s_BV[threadIdx.x].y + Bw[4]);
	ftemp2.x = __frsqrt_rn(taps1.y)*(s_BV[threadIdx.x].x + Bw[1]);     ftemp2.y = __frsqrt_rn(taps1.y)*(s_BV[threadIdx.x].y + Bw[5]);
	//ftemp1.x = __fdividef( (s_BV[threadIdx.x].x + Bw[0]) , c_sqrt_taps[taps1.x]);       ftemp1.y = __fdividef( (s_BV[threadIdx.x].y + Bw[4]) , c_sqrt_taps[taps1.x]);
	//ftemp2.x = __fdividef( (s_BV[threadIdx.x].x + Bw[1]) , c_sqrt_taps[taps1.y]);       ftemp2.y = __fdividef( (s_BV[threadIdx.x].y + Bw[5]) , c_sqrt_taps[taps1.y]);
	if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=startTaps + 2*DIT_value;} else taps1.x=startTaps + DIT_value;
	if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=startTaps + 2*DIT_value;} else taps1.y=startTaps + DIT_value;

	taps2.x=startTaps + 3*DIT_value; taps2.y=startTaps + 4*DIT_value;
	ftemp2.x = __frsqrt_rn(taps2.x)*(s_BV[threadIdx.x].x + Bw[2]);     ftemp2.y = __frsqrt_rn(taps2.x)*(s_BV[threadIdx.x].y + Bw[6]);
	ftemp3.x = __frsqrt_rn(taps2.y)*(s_BV[threadIdx.x].x + Bw[3]);     ftemp3.y = __frsqrt_rn(taps2.y)*(s_BV[threadIdx.x].y + Bw[7]);
	//ftemp2.x = __fdividef( (s_BV[threadIdx.x].x + Bw[2]) , c_sqrt_taps[taps2.x] ); ftemp2.y = __fdividef( (s_BV[threadIdx.x].y + Bw[6]) , c_sqrt_taps[taps2.x] );
	//ftemp3.x = __fdividef( (s_BV[threadIdx.x].x + Bw[3]) , c_sqrt_taps[taps2.y] ); ftemp3.y = __fdividef( (s_BV[threadIdx.x].y + Bw[7]) , c_sqrt_taps[taps2.y] );
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
			ftemp1.x = __frsqrt_rn(taps1.x)*(s_BV[spos].x + Bw[0]);     ftemp1.y = __frsqrt_rn(taps1.x)*(s_BV[spos].y + Bw[4]);
			ftemp2.x = __frsqrt_rn(taps1.y)*(s_BV[spos].x + Bw[1]);     ftemp2.y = __frsqrt_rn(taps1.y)*(s_BV[spos].y + Bw[5]);
			//ftemp1.x = __fdividef( (s_BV[spos].x + Bw[0]) , c_sqrt_taps[taps1.x]);       ftemp1.y = __fdividef( (s_BV[spos].y + Bw[4]) , c_sqrt_taps[taps1.x]);
			//ftemp2.x = __fdividef( (s_BV[spos].x + Bw[1]) , c_sqrt_taps[taps1.y]);       ftemp2.y = __fdividef( (s_BV[spos].y + Bw[5]) , c_sqrt_taps[taps1.y]);
			if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = Taps +2*DIT_value;} else taps1.x = Taps + DIT_value;
			if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = Taps +2*DIT_value;} else taps1.y = Taps + DIT_value;
			
			taps2.x=Taps + 3*DIT_value; taps2.y=Taps + 4*DIT_value;
			ftemp2.x = __frsqrt_rn(taps2.x)*(s_BV[spos].x + Bw[2]);     ftemp2.y = __frsqrt_rn(taps2.x)*(s_BV[spos].y + Bw[6]);
			ftemp3.x = __frsqrt_rn(taps2.y)*(s_BV[spos].x + Bw[3]);     ftemp3.y = __frsqrt_rn(taps2.y)*(s_BV[spos].y + Bw[7]);
			//ftemp2.x = __fdividef( (s_BV[spos].x + Bw[2]) , c_sqrt_taps[taps2.x] ); ftemp2.y = __fdividef( (s_BV[spos].y + Bw[6]) , c_sqrt_taps[taps2.x] );
			//ftemp3.x = __fdividef( (s_BV[spos].x + Bw[3]) , c_sqrt_taps[taps2.y] ); ftemp3.y = __fdividef( (s_BV[spos].y + Bw[7]) , c_sqrt_taps[taps2.y] );
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
	
	if( threadIdx.x<(PD_NTHREADS-(nBoxcars>>1)) ){
		gpos=blockIdx.y*nTimesamples + blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
		d_output_SNR[gpos] = s_SNRs[threadIdx.x].x;
		d_output_SNR[gpos + 1] = s_SNRs[threadIdx.x].y;
		
		d_output_taps[gpos] = s_taps[threadIdx.x].x;
		d_output_taps[gpos + 1] = s_taps[threadIdx.x].y;
		
		d_bv_out[(gpos>>1)] = s_input[threadIdx.x].x;
	}
}


//******************************************************************************************
//******************************************************************************************
//******************************************************************************************


__global__ void PD_GPU_1st_float1_BLN_IF(float const* __restrict__ d_input, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int nTimesamples, int nBoxcars) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, SNR;
	ushort2 taps, taps1, taps2;
	float signal_mean=d_MSD[0];
	float signal_sd=d_MSD[1];
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	
	Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0;
	// Loading data and normalization
	
	gpos = blockIdx.y*nTimesamples + spos;
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

	d_decimated[(gpos>>1)]=Bw[1];
	
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];

	SNR.x=0; SNR.y=0;
	__syncthreads();
	
	ftemp1.x = Bw[0];                    ftemp1.y = Bw[4];
	//ftemp2.x = __frsqrt_rn(2)*Bw[1];     ftemp2.y = __frsqrt_rn(2)*Bw[5];
	ftemp2.x = __fdividef(Bw[1],c_sqrt_taps[2]);     ftemp2.y = __fdividef(Bw[5],c_sqrt_taps[2]);
	
	if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=2;} else taps1.x=1;
	if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=2;} else taps1.y=1;
	
	// Using fast intrinsic
	//ftemp2.x = __frsqrt_rn(3)*Bw[2]; ftemp2.y = __frsqrt_rn(3)*Bw[6];
	//ftemp3.x = __frsqrt_rn(4)*Bw[3]; ftemp3.y = __frsqrt_rn(4)*Bw[7];
	// Using constant memory
	ftemp2.x = __fdividef(Bw[2],c_sqrt_taps[3]);    ftemp2.y = __fdividef(Bw[6],c_sqrt_taps[3]);
	ftemp3.x = __fdividef(Bw[3],c_sqrt_taps[4]);    ftemp3.y = __fdividef(Bw[7],c_sqrt_taps[4]);
	
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
			
			// Using fast intrinsic
			//ftemp1.x = __frsqrt_rn(d+1)*Bw[0];     ftemp1.y = __frsqrt_rn(d+1)*Bw[4];
			//ftemp2.x = __frsqrt_rn(d+2)*Bw[1];     ftemp2.y = __frsqrt_rn(d+2)*Bw[5];
			// Using constant memory
			ftemp1.x = __fdividef(Bw[0],c_sqrt_taps[d+1]);    ftemp1.y = __fdividef(Bw[4],c_sqrt_taps[d+1]);
			ftemp2.x = __fdividef(Bw[1],c_sqrt_taps[d+2]);    ftemp2.y = __fdividef(Bw[5],c_sqrt_taps[d+2]);
			
			if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = d+2;} else taps1.x = d+1;
			if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = d+2;} else taps1.y = d+1;

			// Using fast intrinsic
			//ftemp2.x = __frsqrt_rn(d+3)*Bw[2];     ftemp2.y = __frsqrt_rn(d+3)*Bw[6];
			//ftemp3.x = __frsqrt_rn(d+4)*Bw[3];     ftemp3.y = __frsqrt_rn(d+4)*Bw[7];
			// Using constant memory
			ftemp2.x = __fdividef(Bw[2],c_sqrt_taps[d+3]);    ftemp2.y = __fdividef(Bw[6],c_sqrt_taps[d+3]);
			ftemp3.x = __fdividef(Bw[3],c_sqrt_taps[d+4]);    ftemp3.y = __fdividef(Bw[7],c_sqrt_taps[d+4]);
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
		
		d_bv_out[(gpos>>1)] = s_input[threadIdx.x].x;
	}
}


__global__ void PD_GPU_Nth_float1_BLN_IF(float const* __restrict__ d_input, float *d_bv_in, float *d_bv_out, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int nTimesamples, int nBoxcars, int startTaps, int DIT_value, int shift) {
	__shared__ float2 s_input[PD_NTHREADS];
	__shared__ float2 s_BV[PD_NTHREADS];
	__shared__ float2 s_SNRs[PD_NTHREADS];
	__shared__ ushort2 s_taps[PD_NTHREADS];
	
	
	int d, gpos, spos;
	float Bw[8];
	float2 ftemp1, ftemp2, ftemp3, SNR;
	ushort2 taps, taps1, taps2;
	int Taps;
	
	Bw[0]=0; Bw[1]=0; Bw[2]=0; Bw[3]=0; Bw[4]=0; Bw[5]=0; Bw[6]=0; Bw[7]=0;
	s_BV[threadIdx.x].x=0; s_BV[threadIdx.x].y=0;
	
	spos=blockIdx.x*(2*PD_NTHREADS-nBoxcars) + 2*threadIdx.x;
	
	gpos=blockIdx.y*nTimesamples + spos;
	ftemp1.x=d_input[gpos+shift];
	ftemp1.y=d_input[gpos+shift + 1];
	ftemp2.x=d_input[gpos+shift + 2];
	ftemp2.y=d_input[gpos+shift + 3];
	ftemp3.x=d_input[gpos+shift + 4];
	
	
	Bw[0]=ftemp1.x; Bw[4]=ftemp1.y;
	Bw[1]=ftemp1.x+ftemp1.y; Bw[5]=ftemp1.y+ftemp2.x;
	Bw[2]=Bw[1] + ftemp2.x; Bw[6]=Bw[5] + ftemp2.y;
	Bw[3]=Bw[1] + ftemp2.x + ftemp2.y; Bw[7] = Bw[5] + ftemp2.y + ftemp3.x;
	
	
	d_decimated[(gpos>>1)]=Bw[1];
    
	
	s_BV[threadIdx.x].x=d_bv_in[gpos];	
	s_BV[threadIdx.x].y=d_bv_in[gpos+1];
	
	
	s_input[threadIdx.x].x = Bw[3];
	s_input[threadIdx.x].y = Bw[7];
	
	SNR.x=0; SNR.y=0;
	__syncthreads();
	
	taps1.x = startTaps + DIT_value; taps1.y=startTaps + 2*DIT_value;
	ftemp1.x = __frsqrt_rn(taps1.x)*(s_BV[threadIdx.x].x + Bw[0]);     ftemp1.y = __frsqrt_rn(taps1.x)*(s_BV[threadIdx.x].y + Bw[4]);
	ftemp2.x = __frsqrt_rn(taps1.y)*(s_BV[threadIdx.x].x + Bw[1]);     ftemp2.y = __frsqrt_rn(taps1.y)*(s_BV[threadIdx.x].y + Bw[5]);
	//ftemp1.x = __fdividef( (s_BV[threadIdx.x].x + Bw[0]) , c_sqrt_taps[taps1.x]);       ftemp1.y = __fdividef( (s_BV[threadIdx.x].y + Bw[4]) , c_sqrt_taps[taps1.x]);
	//ftemp2.x = __fdividef( (s_BV[threadIdx.x].x + Bw[1]) , c_sqrt_taps[taps1.y]);       ftemp2.y = __fdividef( (s_BV[threadIdx.x].y + Bw[5]) , c_sqrt_taps[taps1.y]);
	if(ftemp2.x>ftemp1.x) {ftemp1.x=ftemp2.x; taps1.x=startTaps + 2*DIT_value;} else taps1.x=startTaps + DIT_value;
	if(ftemp2.y>ftemp1.y) {ftemp1.y=ftemp2.y; taps1.y=startTaps + 2*DIT_value;} else taps1.y=startTaps + DIT_value;

	taps2.x=startTaps + 3*DIT_value; taps2.y=startTaps + 4*DIT_value;
	ftemp2.x = __frsqrt_rn(taps2.x)*(s_BV[threadIdx.x].x + Bw[2]);     ftemp2.y = __frsqrt_rn(taps2.x)*(s_BV[threadIdx.x].y + Bw[6]);
	ftemp3.x = __frsqrt_rn(taps2.y)*(s_BV[threadIdx.x].x + Bw[3]);     ftemp3.y = __frsqrt_rn(taps2.y)*(s_BV[threadIdx.x].y + Bw[7]);
	//ftemp2.x = __fdividef( (s_BV[threadIdx.x].x + Bw[2]) , c_sqrt_taps[taps2.x] ); ftemp2.y = __fdividef( (s_BV[threadIdx.x].y + Bw[6]) , c_sqrt_taps[taps2.x] );
	//ftemp3.x = __fdividef( (s_BV[threadIdx.x].x + Bw[3]) , c_sqrt_taps[taps2.y] ); ftemp3.y = __fdividef( (s_BV[threadIdx.x].y + Bw[7]) , c_sqrt_taps[taps2.y] );
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
			ftemp1.x = __frsqrt_rn(taps1.x)*(s_BV[spos].x + Bw[0]);     ftemp1.y = __frsqrt_rn(taps1.x)*(s_BV[spos].y + Bw[4]);
			ftemp2.x = __frsqrt_rn(taps1.y)*(s_BV[spos].x + Bw[1]);     ftemp2.y = __frsqrt_rn(taps1.y)*(s_BV[spos].y + Bw[5]);
			//ftemp1.x = __fdividef( (s_BV[spos].x + Bw[0]) , c_sqrt_taps[taps1.x]);       ftemp1.y = __fdividef( (s_BV[spos].y + Bw[4]) , c_sqrt_taps[taps1.x]);
			//ftemp2.x = __fdividef( (s_BV[spos].x + Bw[1]) , c_sqrt_taps[taps1.y]);       ftemp2.y = __fdividef( (s_BV[spos].y + Bw[5]) , c_sqrt_taps[taps1.y]);
			if(ftemp2.x>ftemp1.x) {ftemp1.x = ftemp2.x; taps1.x = Taps +2*DIT_value;} else taps1.x = Taps + DIT_value;
			if(ftemp2.y>ftemp1.y) {ftemp1.y = ftemp2.y; taps1.y = Taps +2*DIT_value;} else taps1.y = Taps + DIT_value;
			
			taps2.x=Taps + 3*DIT_value; taps2.y=Taps + 4*DIT_value;
			ftemp2.x = __frsqrt_rn(taps2.x)*(s_BV[spos].x + Bw[2]);     ftemp2.y = __frsqrt_rn(taps2.x)*(s_BV[spos].y + Bw[6]);
			ftemp3.x = __frsqrt_rn(taps2.y)*(s_BV[spos].x + Bw[3]);     ftemp3.y = __frsqrt_rn(taps2.y)*(s_BV[spos].y + Bw[7]);
			//ftemp2.x = __fdividef( (s_BV[spos].x + Bw[2]) , c_sqrt_taps[taps2.x] ); ftemp2.y = __fdividef( (s_BV[spos].y + Bw[6]) , c_sqrt_taps[taps2.x] );
			//ftemp3.x = __fdividef( (s_BV[spos].x + Bw[3]) , c_sqrt_taps[taps2.y] ); ftemp3.y = __fdividef( (s_BV[spos].y + Bw[7]) , c_sqrt_taps[taps2.y] );
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
		
		d_bv_out[(gpos>>1)] = s_input[threadIdx.x].x;
	}
}



#endif
