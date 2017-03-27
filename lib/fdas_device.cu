/* Device functions for acceleration search  */
/* 06/11/2014 Sofia Dimoudi sofia.dimoudi@oerc.ox.ac.uk */


#include <stdio.h>
#include <cuda_runtime.h>
#include "headers/fdas_device.h"
#include "headers/params.h"

static __device__ __inline__ float2 Get_W_value(int N, int m){
	float2 ctemp;
	ctemp.x=cosf( -2.0f*3.141592654f*fdividef( (float) m, (float) N) );
	ctemp.y=sinf( -2.0f*3.141592654f*fdividef( (float) m, (float) N) );	
	return(ctemp);
}

static __device__ __inline__ float2 Get_W_value_inverse(int N, int m){
	float2 ctemp;
	ctemp.x=cosf(2.0f*3.141592654f*fdividef( (float)(m), (float)(N)) );
	ctemp.y=sinf( 2.0f*3.141592654f*fdividef( (float)(m), (float)(N)) );
	return(ctemp);
}
/*
static __device__ __inline__ float2 Get_W_value_float(float N, float m){
	float2 ctemp;
	ctemp.x=-cosf( 6.283185f*fdividef( m, N) - 3.141592654f );
	ctemp.y=sinf( 6.283185f*fdividef( m, N) - 3.141592654f );
	return(ctemp);
}
*/
static __device__ __inline__ void do_FFT_no_reorder(float2 *s_input, int N, int bits){ // in-place
	float2 A_DFT_value, B_DFT_value, ftemp2, ftemp;
	float2 WB;
	
	int r, Bj, Bk, PoTm1, A_read_index, B_read_index, A_write_index, B_write_index, Nhalf;

	Nhalf=N>>1;
	
	//-----> FFT
	//--> 
	PoTm1=1;
	
	A_read_index=threadIdx.x;
	B_read_index=threadIdx.x + Nhalf;
		
	A_write_index=2*threadIdx.x;
	B_write_index=2*threadIdx.x+1;
	
	for(r=1;r<bits;r++){
		if(threadIdx.x<Nhalf){
			Bj=(2*threadIdx.x+1)>>r;

			Bk=PoTm1*Bj;
			
			// A element
			ftemp  = s_input[A_read_index];
			ftemp2 = s_input[B_read_index];
			
			// WA=1
			
			A_DFT_value.x=ftemp.x + ftemp2.x;
			A_DFT_value.y=ftemp.y + ftemp2.y;
			
			WB = Get_W_value(N,Bk);
			
			B_DFT_value.x=WB.x*(ftemp.x - ftemp2.x) - WB.y*(ftemp.y - ftemp2.y);
			B_DFT_value.y=WB.x*(ftemp.y - ftemp2.y) + WB.y*(ftemp.x - ftemp2.x);
			
			PoTm1=PoTm1<<1;
		}
		__syncthreads();
		if(threadIdx.x<Nhalf){
			s_input[A_write_index]=A_DFT_value;
			s_input[B_write_index]=B_DFT_value;
		}
		__syncthreads();
	}
	
	//-----------------------------------------------
	//----- Last exchange
	if(threadIdx.x<Nhalf){
		ftemp  = s_input[A_read_index];
		ftemp2 = s_input[B_read_index];
		
		A_DFT_value.x = ftemp.x + ftemp2.x;
		A_DFT_value.y = ftemp.y + ftemp2.y;
		B_DFT_value.x = ftemp.x - ftemp2.x;
		B_DFT_value.y = ftemp.y - ftemp2.y;
	}
	__syncthreads();
	if(threadIdx.x<Nhalf){
		s_input[A_write_index]=A_DFT_value;
		s_input[B_write_index]=B_DFT_value;
	}
	__syncthreads();
	//----- Last exchange
	//-----------------------------------------------
		
	//-------> END
}

static __device__ __inline__ void do_IFFT_no_reorder(float2 *s_input, int N, int bits){ // in-place
	float2 A_DFT_value, B_DFT_value;
	float2 W;
	float2 ftemp, ftemp2;

	int local_id, warp_id;
	int j, m_param;
	int parity, itemp;
	int A_read_index,B_read_index;
	int PoT, PoTp1, q;
	int Nhalf;

	Nhalf=N>>1;
	
	local_id = threadIdx.x & (WARP - 1);
	warp_id = threadIdx.x/WARP;
	
	//-----> FFT
	//-->
	if(threadIdx.x<Nhalf){
		PoT=1;
		PoTp1=2;	

		//--> First iteration
		itemp=local_id - (local_id&4294967294);
		parity=(1-itemp*2);
		A_DFT_value=s_input[local_id + warp_id*2*WARP];
		B_DFT_value=s_input[local_id + warp_id*2*WARP + WARP];
		
		A_DFT_value.x=parity*A_DFT_value.x+ __shfl(A_DFT_value.x,local_id+parity);
		A_DFT_value.y=parity*A_DFT_value.y+ __shfl(A_DFT_value.y,local_id+parity);
		
		B_DFT_value.x=parity*B_DFT_value.x+ __shfl(B_DFT_value.x,local_id+parity);
		B_DFT_value.y=parity*B_DFT_value.y+ __shfl(B_DFT_value.y,local_id+parity);
		
		//--> First iteration
		
		PoT=2;
		PoTp1=4;
		
		for(q=2;q<6;q++){
			
			m_param = local_id & (PoTp1 - 1);
			itemp=m_param/PoT;
			parity=(1-itemp*2);
			
			W=Get_W_value_inverse(PoT*2,m_param);

			A_read_index=local_id+parity*itemp*PoT;
			B_read_index=local_id+(1-itemp)*PoT;
			
			ftemp2.x=__shfl(A_DFT_value.x,B_read_index);
			ftemp2.y=__shfl(A_DFT_value.y,B_read_index);					
			A_DFT_value.x=__shfl(A_DFT_value.x,A_read_index) + W.x*ftemp2.x-W.y*ftemp2.y;
			A_DFT_value.y=__shfl(A_DFT_value.y,A_read_index) + W.x*ftemp2.y+W.y*ftemp2.x;
			
			ftemp.x=__shfl(B_DFT_value.x,B_read_index);
			ftemp.y=__shfl(B_DFT_value.y,B_read_index);					
			B_DFT_value.x=__shfl(B_DFT_value.x,A_read_index) + W.x*ftemp.x - W.y*ftemp.y;
			B_DFT_value.y=__shfl(B_DFT_value.y,A_read_index) + W.x*ftemp.y + W.y*ftemp.x;		
			
			PoT=PoT<<1;
			PoTp1=PoTp1<<1;
		}	
		
		s_input[local_id + warp_id*2*WARP]=A_DFT_value;
		s_input[local_id + warp_id*2*WARP + WARP]=B_DFT_value;
	}
	__syncthreads();
	
	for(q=6;q<=bits;q++){
		if(threadIdx.x<Nhalf){
			m_param = threadIdx.x & (PoT - 1);
			j=threadIdx.x>>(q-1);
			
			W=Get_W_value_inverse(PoTp1,m_param);

			A_read_index=j*PoTp1 + m_param;
			B_read_index=j*PoTp1 + m_param + PoT;
			
			ftemp  = s_input[A_read_index];
			ftemp2 = s_input[B_read_index];
			
			A_DFT_value.x=ftemp.x + W.x*ftemp2.x - W.y*ftemp2.y;
			A_DFT_value.y=ftemp.y + W.x*ftemp2.y + W.y*ftemp2.x;
			
			B_DFT_value.x=ftemp.x - W.x*ftemp2.x + W.y*ftemp2.y;
			B_DFT_value.y=ftemp.y - W.x*ftemp2.y - W.y*ftemp2.x;

			PoT=PoT<<1;
			PoTp1=PoTp1<<1;		
		}
		__syncthreads();
		if(threadIdx.x<Nhalf){
			s_input[A_read_index]=A_DFT_value;
			s_input[B_read_index]=B_DFT_value;
		}
		__syncthreads();
		

	}
}


/* ----------------------------------- */

static __device__ inline float pwcalc(float2 cpx)
{
  float pw = (cpx.x*cpx.x + cpx.y*cpx.y);
  return pw;

}

/* -------- FDAS KERNELS ----------*/

__global__ void cuda_overlap_copy(float2* d_ext_data, float2* d_cpx_signal, int sigblock, int sig_rfftlen, int sig_tot_convlen, int kern_offset, int total_blocks)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int numthreads = blockDim.x * gridDim.x; 

  //initialize the array  
  for (int i = tid; i < sig_tot_convlen; i += numthreads){
    d_ext_data[i].x = 0.0f;
    d_ext_data[i].y =0.0f;
  }

  if (tid >= kern_offset && tid < KERNLEN ) //firdt block
    d_ext_data[tid] = d_cpx_signal[tid - kern_offset ];
  
  for (int i = 1; i < total_blocks; i++ ){ // copy overlapped blocks
    d_ext_data[( i*KERNLEN + tid) ] = d_cpx_signal[i*sigblock - kern_offset + tid ];
  } 

}

__global__ void cuda_overlap_copy_smallblk(float2* d_ext_data, float2* d_cpx_signal, int sigblock, int sig_rfftlen, int sig_tot_convlen, int kern_offset, int total_blocks)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int read_idx = blockIdx.x*sigblock - kern_offset + threadIdx.x;
  int write_idx = blockIdx.x*KERNLEN + threadIdx.x;

  //initialize the array  
  if (tid < sig_tot_convlen){
    d_ext_data[tid].x = 0.0f;
    d_ext_data[tid].y =0.0f;
  }

  if (threadIdx.x >= kern_offset && blockIdx.x == 0 ) //first block
    d_ext_data[threadIdx.x] = d_cpx_signal[threadIdx.x - kern_offset ];
  
  // copy overlapped blocks
  if (blockIdx.x > 0 && read_idx < sig_rfftlen){    
    d_ext_data[write_idx] = d_cpx_signal[read_idx];
  }
}

__global__ void cuda_convolve_reg_1d_halftemps(float2* d_kernel, float2* d_signal, float2* d_ffdot_plane,int sig_tot_convlen, float scale)
{
  /* Load only half templates - use register arrays */
  int tx = threadIdx.x;
  int bx = blockIdx.x;
  int tidx = bx* blockDim.x + tx;

  float2 tempkern[ZMAX/2]; // half template array, 1 column per thread
  for (int i=0; i < ZMAX/2; i++ )
    tempkern[i] = __ldg(&d_kernel[i*KERNLEN + tidx]) ;

  for (int j = 0; j < sig_tot_convlen; j+=KERNLEN){
    float2 local_data = __ldg(&d_signal[tidx + j]);
#pragma unroll
    for (int i = 0; i < ZMAX/2; i++){
      //upper half
      d_ffdot_plane[i*sig_tot_convlen + j + tidx ].x = (local_data.x*tempkern[i].x - local_data.y*tempkern[i].y) * scale;
      d_ffdot_plane[i*sig_tot_convlen + j + tidx ].y = (local_data.x*tempkern[i].y + local_data.y*tempkern[i].x) * scale;
      //complex conjugate filter lower half
      d_ffdot_plane[(ZMAX - i)*sig_tot_convlen + j + tidx ].x = (local_data.x*tempkern[i].x + local_data.y*tempkern[i].y) * scale;
      d_ffdot_plane[(ZMAX - i)*sig_tot_convlen + j + tidx ].y = (-local_data.x*tempkern[i].y + local_data.y*tempkern[i].x) * scale;
    }

    //do z=0
    float2 tempz = __ldg(&d_kernel[(ZMAX/2)*KERNLEN + tidx]) ;
    d_ffdot_plane[(ZMAX/2)*sig_tot_convlen + j + tidx ].x = (local_data.x*tempz.x - local_data.y*tempz.y) * scale;
    d_ffdot_plane[(ZMAX/2)*sig_tot_convlen + j + tidx ].y = (local_data.x*tempz.y + local_data.y*tempz.x) * scale;
  }
}

__global__ void cuda_ffdotpow_concat_2d(float2* d_ffdot_plane_cpx, float* d_ffdot_plane, int sigblock, int kern_offset, int total_blocks, int sig_tot_convlen, int sig_totlen)
 /* Copies useful points from convolved complex f-fdot array (discarding contaminated ends) */
 /* calculates complex signal power and places result in float f-fdot array */
 {
   int tx = threadIdx.x;
   int ty = threadIdx.y;
   int bx = blockIdx.x;
   int by = blockIdx.y;
   int tidx = bx* blockDim.x + tx;
   int tidy = by* blockDim.y + ty;
     
   for (int i = 0; i < total_blocks; ++i){
     if (tidx < sigblock){
       float2  local_data =  __ldg(&d_ffdot_plane_cpx[(tidy * sig_tot_convlen) + i*KERNLEN  + tidx + kern_offset ]) ;
       d_ffdot_plane[ (tidy * sig_totlen) + i*sigblock + tidx] = (local_data.x*local_data.x + local_data.y*local_data.y); 
     }
   }
 }

 __global__ void cuda_ffdotpow_concat_2d_inbin(float2* d_ffdot_plane_cpx, float* d_ffdot_plane, int sigblock, int kern_offset, int total_blocks, int sig_tot_convlen, int sig_totlen)
 /* Copies useful points from convolved complex f-fdot array (discarding contaminated ends) */
 /* calculates complex signal power and places result in float f-fdot array */
 /* This one performs interbinning on the complex signal before power calculations */
 {
   int tx = threadIdx.x;
   int bx = blockIdx.x;
   int by = blockIdx.y;
   int tidx = bx* blockDim.x + tx;

   __shared__ float2 local_data[PTBSIZEX + 1];
   __shared__ float local_data_inbin[2*PTBSIZEX];
 
   for (int i = 0; i < total_blocks; ++i){
     //     __syncthreads();
     if (tidx < sigblock){
	 //first load complex f-fdot blocks between overlap offsets
       // onto shared memory
       local_data[tx] =  __ldg(&d_ffdot_plane_cpx[(by * sig_tot_convlen) + i*KERNLEN  + tidx + kern_offset  ]) ;
       __syncthreads();

       //load one extra point in the x direction
       if(tx==PTBSIZEX-1)
	 local_data[tx+1] = __ldg(&d_ffdot_plane_cpx[ (by * sig_tot_convlen) + i*KERNLEN  + tidx + kern_offset + 1 ]) ;
       __syncthreads();

       // pick first points from the next segment, in i+1
       if (tidx == sigblock - 1){
	 //	 if(tx==PTBSIZEX-1)
	   local_data[tx+1] = __ldg(&d_ffdot_plane_cpx[(by * sig_tot_convlen) + (i+1)*KERNLEN  + kern_offset ]) ;
       }
       __syncthreads();

	 //spread signal in interleaving bins within each block, and write power to each even bin
       local_data_inbin[2*tx] = local_data[tx].x*local_data[tx].x + local_data[tx].y*local_data[tx].y; 
       __syncthreads();

	 // fill in the empty bins: 
       // P_k+1/2 = |F_k+1/2|^2 =~|pi/4 * (F_k - F_k-1)|^2 
       // (see Handbook of Pulsar astronomy p.135)
       local_data_inbin[2*tx + 1] =  0.616850275f * ((local_data[tx].x - local_data[tx+1].x) * (local_data[tx].x - local_data[tx+1].x) + (local_data[tx].y - local_data[tx+1].y) * (local_data[tx].y - local_data[tx+1].y));
       __syncthreads();
     
       //write data back to global memory contiguously
       unsigned int inbin_idx = (unsigned int)(2*by*sig_totlen + i*2*sigblock + 2*bx*blockDim.x) + tx;
       d_ffdot_plane[inbin_idx] = local_data_inbin[tx]; 
       d_ffdot_plane[inbin_idx + PTBSIZEX] = local_data_inbin[tx + PTBSIZEX]; 
     }
   }
 }

__global__ void cuda_ffdotpow_concat_2d_ndm_inbin(float2* d_ffdot_plane_cpx, float* d_ffdot_plane, int kernlen, int siglen, int nkern, int kern_offset, int total_blocks, int sig_tot_convlen, int sig_totlen, int ndm)
/* Copies useful points from convolved complex f-fdot array (discarding contaminated ends) */
/* calculates complex signal power and places result in float f-fdot array */
/* This one performs interbinning on the complex signal before power calculations */
{
  int tx = threadIdx.x;
  int bx = blockIdx.x;
  int by = blockIdx.y;
  int tidx = bx* blockDim.x + tx;
  //int xthreads = blockDim.x * gridDim.x;
  unsigned int ffdotlen = NKERN * sig_totlen;
  unsigned int ffdotlencpx = NKERN * sig_tot_convlen;
  //   int inbin = 2; //temporary, for test

  
  __shared__ float2 local_data[PTBSIZEX + 1];
  __shared__ float local_data_inbin[2*PTBSIZEX];
 
  for (int dm = 0; dm < ndm; dm++){
    unsigned int dmidx_cpx = dm*ffdotlencpx;
    unsigned int dmidx_f = dm*ffdotlen;
    for (int i = 0; i < total_blocks; i++){
      //       __syncthreads();
      if (tidx < siglen){
	//first load complex f-fdot blocks between overlap offsets onto shared memory
	local_data[tx] =  __ldg(&d_ffdot_plane_cpx[dmidx_cpx + (by * sig_tot_convlen) + i*kernlen  + tidx + kern_offset ]) ;
	__syncthreads();
	//load one extra point in the x direction
	if(tx==PTBSIZEX-1)
	    local_data[tx+1] = __ldg(&d_ffdot_plane_cpx[dmidx_cpx + (by * sig_tot_convlen) + i*kernlen  + tidx + kern_offset + 1 ]) ;
	__syncthreads();
	// pick first points from the next segment, in i+1
	if (tidx == siglen - 1){ 
	  // if(tx==PTBSIZEX-1)
	   local_data[tx+1] = __ldg(&d_ffdot_plane_cpx[dmidx_cpx + (by * sig_tot_convlen) + (i+1)*kernlen  + kern_offset ]) ;
	}
	__syncthreads();
	//spread signal in interleaving bins within each block, and write power to each even bin
	local_data_inbin[2*tx] = local_data[tx].x*local_data[tx].x + local_data[tx].y*local_data[tx].y; 
	__syncthreads();
	// fill in the empty bins: P_k+1/2 = |F_k+1/2|^2 =~|pi/4 * (F_k - F_k-1)|^2 (see Handbook of Pulsar astronomy p.135)
       
	local_data_inbin[2*tx + 1] =  0.616850275f * ((local_data[tx].x - local_data[tx+1].x) * (local_data[tx].x - local_data[tx+1].x) + (local_data[tx].y - local_data[tx+1].y) * (local_data[tx].y - local_data[tx+1].y));
	__syncthreads();
	//write data back to global memory
	//	int inbin_idx = 2*dmidx_f + 2*(by * sig_totlen) + i*2*siglen + 2*bx* blockDim.x + tx ; 
	//	if (inbin_idx + PTBSIZEX < ffdotlen){
	  d_ffdot_plane[2*dmidx_f + 2*(by * sig_totlen) + i*2*siglen + 2*bx* blockDim.x + tx] = local_data_inbin[tx];
	  d_ffdot_plane[ 2*dmidx_f + 2*(by * sig_totlen) + i*2*siglen + 2*bx* blockDim.x + tx + PTBSIZEX] = local_data_inbin[tx+PTBSIZEX]; 
	  //	}
	// __syncthreads();
      }
    }
  }
}

#ifndef NOCUST
__global__ void customfft_fwd_temps_no_reorder(float2* d_signal)
{
  int tx = threadIdx.x;
  int bx = blockIdx.x;

  __shared__ float2 s_input[KERNLEN]; 

 s_input[tx]=d_signal[tx + bx*KERNLEN];
  __syncthreads();

  //call custom device fft
  do_FFT_no_reorder(s_input,KERNLEN,NEXP);
  __syncthreads();

  d_signal[tx +bx*KERNLEN] = s_input[tx];

}

__global__ void cuda_convolve_customfft_wes_no_reorder02(float2* d_kernel, float2* d_signal, float *d_ffdot_pw, int sigblock, int sig_tot_convlen, int sig_totlen, int offset, float scale)
/* convolution kernel using Karel Adamek's custom FFT, deployed here with modifications by Wes Armour.
   It performs the forward FFT and then loops through filters.(1-d blocks)
   It also uses half the templates and computes the rest using the complex conjugate.
   The modifications are optimizing speed by computing FFTs on two templates with one synchronization point. 
   In this version we apply Karel's kernel with removed de-shuffling of the data during the FFT.*/
{
  int tx = threadIdx.x;
  int bx = blockIdx.x;
  float2 ffdotcpx;
  float2 ffdotcpx_trans;
  float2 local_data; //fft'd data

  float one;
  float two;
  float three;
  float four;

  float2 A_DFT_value, B_DFT_value, sA_DFT_value, sB_DFT_value;
  float2 ftemp2, ftemp, stemp2, stemp;
  float2 W;
  int local_id, warp_id;
  int j, m_param;
  //int load_id, i, n;
  int parity, itemp;
  int A_read_index,B_read_index;
  //  int A_write_index,B_write_index;
  int PoT, PoTp1, q;
	
  __shared__ float2 s_input[KERNLEN]; //signal input data to FFT
  __shared__ float2 s_input_trans[KERNLEN]; //static allocation


  int index =  bx*sigblock + tx;
  int tidx =  bx * blockDim.x + tx;

  //load signal data to shared memory
  s_input[tx]=__ldg(&d_signal[tidx]);
  __syncthreads();


  //call custom device fft
  do_FFT_no_reorder(s_input,KERNLEN,NEXP);

  local_data.x = s_input[tx].x;
  local_data.y = s_input[tx].y;
  // __syncthreads();

  //complex multiplication power calculation loop over template columns
  #pragma unroll 
  for (int i = 0; i < ZMAX/2; i++){
    //__syncthreads();
    // fourier domain convolution
    one = local_data.x*__ldg(&d_kernel[i*KERNLEN + tx].x);
    two = local_data.y*__ldg(&d_kernel[i*KERNLEN + tx].y);
    three = local_data.x*__ldg(&d_kernel[i*KERNLEN + tx].y);
    four = local_data.y*__ldg(&d_kernel[i*KERNLEN + tx].x);

    ffdotcpx.x = ((one - two) * scale);
    ffdotcpx.y = ((three + four) * scale);
    s_input[tx] = ffdotcpx;

    ffdotcpx_trans.x = ((one + two) * scale);
    ffdotcpx_trans.y = ((four -three) * scale);
    s_input_trans[tx] = ffdotcpx_trans;
    __syncthreads();
    //---------

    //inverse fft
    local_id = threadIdx.x & (WARP - 1);
    warp_id = threadIdx.x/WARP;
    //-----> FFT
    //-->
    if(threadIdx.x<KERNLEN/2){
      PoT=1;
      PoTp1=2;	
      //--> First iteration
      itemp=local_id - (local_id&4294967294);
      parity=(1-itemp*2);
      //1st template
      A_DFT_value=s_input[local_id + warp_id*2*WARP];
      sA_DFT_value=s_input_trans[local_id + warp_id*2*WARP];
      B_DFT_value=s_input[local_id + warp_id*2*WARP + WARP];
      sB_DFT_value=s_input_trans[local_id + warp_id*2*WARP + WARP];

      //1st template
      A_DFT_value.x=parity*A_DFT_value.x+ __shfl(A_DFT_value.x,local_id+parity);
      A_DFT_value.y=parity*A_DFT_value.y+ __shfl(A_DFT_value.y,local_id+parity);
      B_DFT_value.x=parity*B_DFT_value.x+ __shfl(B_DFT_value.x,local_id+parity);
      B_DFT_value.y=parity*B_DFT_value.y+ __shfl(B_DFT_value.y,local_id+parity);
      //2nd template

      sA_DFT_value.x=parity*sA_DFT_value.x+ __shfl(sA_DFT_value.x,local_id+parity);
      sA_DFT_value.y=parity*sA_DFT_value.y+ __shfl(sA_DFT_value.y,local_id+parity);
      sB_DFT_value.x=parity*sB_DFT_value.x+ __shfl(sB_DFT_value.x,local_id+parity);
      sB_DFT_value.y=parity*sB_DFT_value.y+ __shfl(sB_DFT_value.y,local_id+parity);

      //--> First iteration end
		
      PoT=2;
      PoTp1=4;
      float fPoT=2.0f;

      for(q=2;q<6;q++){
	m_param = local_id & (PoTp1 - 1);
	//itemp=m_param/PoT;
	itemp=__float2int_rz(__fdividef(__int2float_rz(m_param),fPoT));
	parity=(1-itemp*2);
	  
	//	W=Get_W_value_inverse(PoT*2,m_param);
	W.x=cosf( 2.0f*3.141592654f*fdividef( (float) m_param, (float) (PoT*2)) );
	W.y=sinf( 2.0f*3.141592654f*fdividef( (float) m_param, (float) (PoT*2)) );
	  
	A_read_index=local_id+parity*itemp*PoT;
	B_read_index=local_id+(1-itemp)*PoT;
			
	//1st template
	ftemp2.x=__shfl(A_DFT_value.x,B_read_index);
	ftemp2.y=__shfl(A_DFT_value.y,B_read_index);				
	//2nd template
	stemp2.x=__shfl(sA_DFT_value.x,B_read_index);
	stemp2.y=__shfl(sA_DFT_value.y,B_read_index);
	
	//1st template
	A_DFT_value.x=__shfl(A_DFT_value.x,A_read_index) + W.x*ftemp2.x-W.y*ftemp2.y;
	A_DFT_value.y=__shfl(A_DFT_value.y,A_read_index) + W.x*ftemp2.y+W.y*ftemp2.x;
	//2nd template
	sA_DFT_value.x=__shfl(sA_DFT_value.x,A_read_index) + W.x*stemp2.x-W.y*stemp2.y;
	sA_DFT_value.y=__shfl(sA_DFT_value.y,A_read_index) + W.x*stemp2.y+W.y*stemp2.x;
	

	//1st template	  
	ftemp.x=__shfl(B_DFT_value.x,B_read_index);
	ftemp.y=__shfl(B_DFT_value.y,B_read_index);
	//2nd template	  
	stemp.x=__shfl(sB_DFT_value.x,B_read_index);
	stemp.y=__shfl(sB_DFT_value.y,B_read_index);
	
	//1st template
	B_DFT_value.x=__shfl(B_DFT_value.x,A_read_index) + W.x*ftemp.x - W.y*ftemp.y;
	B_DFT_value.y=__shfl(B_DFT_value.y,A_read_index) + W.x*ftemp.y + W.y*ftemp.x;		
	//2nd template
	sB_DFT_value.x=__shfl(sB_DFT_value.x,A_read_index) + W.x*stemp.x - W.y*stemp.y;
	sB_DFT_value.y=__shfl(sB_DFT_value.y,A_read_index) + W.x*stemp.y + W.y*stemp.x;		
	  
	PoT=PoT<<1;
	fPoT=fPoT*2.0f;
	PoTp1=PoTp1<<1;
      }
      //1st template
      s_input[local_id + warp_id*2*WARP]=A_DFT_value;
      s_input[local_id + warp_id*2*WARP + WARP]=B_DFT_value;
      //2nd template
      s_input_trans[local_id + warp_id*2*WARP]=sA_DFT_value;
      s_input_trans[local_id + warp_id*2*WARP + WARP]=sB_DFT_value;
      
    }
//    __syncthreads();
      
    for(q=6;q<=NEXP;q++){
      if(threadIdx.x<KERNLEN/2){
	m_param = threadIdx.x & (PoT - 1);
	j=threadIdx.x>>(q-1);
	  
	W=Get_W_value_inverse(PoTp1,m_param);

	A_read_index=j*PoTp1 + m_param;
	B_read_index=j*PoTp1 + m_param + PoT;
	  
	//1st template
	ftemp  = s_input[A_read_index];
	ftemp2 = s_input[B_read_index];
	//2nd template
	stemp  = s_input_trans[A_read_index];
	stemp2 = s_input_trans[B_read_index];
	
	//1st template
	A_DFT_value.x=ftemp.x + W.x*ftemp2.x - W.y*ftemp2.y;
	A_DFT_value.y=ftemp.y + W.x*ftemp2.y + W.y*ftemp2.x;
	//2nd template
	sA_DFT_value.x=stemp.x + W.x*stemp2.x - W.y*stemp2.y;
	sA_DFT_value.y=stemp.y + W.x*stemp2.y + W.y*stemp2.x;
	
	//1st template	  
	B_DFT_value.x=ftemp.x - W.x*ftemp2.x + W.y*ftemp2.y;
	B_DFT_value.y=ftemp.y - W.x*ftemp2.y - W.y*ftemp2.x;
	//2nd template
	sB_DFT_value.x=stemp.x - W.x*stemp2.x + W.y*stemp2.y;
	sB_DFT_value.y=stemp.y - W.x*stemp2.y - W.y*stemp2.x;
	
	PoT=PoT<<1;
	PoTp1=PoTp1<<1;		
     // }
      // __syncthreads();
     // if(threadIdx.x<KERNLEN/2){
	//1st template	  
	s_input[A_read_index]=A_DFT_value;
	s_input[B_read_index]=B_DFT_value;
	//2nd template	  
	s_input_trans[A_read_index]=sA_DFT_value;
	s_input_trans[B_read_index]=sB_DFT_value;
	
      }
    }
        
       __syncthreads();
    if(tx < sigblock){
      d_ffdot_pw[i*sig_totlen + index] = pwcalc(s_input[tx +offset]);

      d_ffdot_pw[(ZMAX - i )*sig_totlen + index] = pwcalc(s_input_trans[tx+offset]);
    }
  }

  // now do z=0
  one = local_data.x*__ldg(&d_kernel[(ZMAX/2)*KERNLEN + tx].x);
  two = local_data.y*__ldg(&d_kernel[(ZMAX/2)*KERNLEN + tx].y);
  three = local_data.x*__ldg(&d_kernel[(ZMAX/2)*KERNLEN + tx].y);
  four = local_data.y*__ldg(&d_kernel[(ZMAX/2)*KERNLEN + tx].x);
  
  ffdotcpx.x = ((one - two) * scale);
  ffdotcpx.y = ((three + four) * scale);
  s_input[tx] = ffdotcpx;
  __syncthreads();
  //call inverse fft
  do_IFFT_no_reorder(s_input, KERNLEN ,NEXP);

  // write powers 
  if(tx < sigblock)
    d_ffdot_pw[(ZMAX/2) * sig_totlen + index] = pwcalc(s_input[tx+offset]);

  //-------END
}

__global__ void cuda_convolve_customfft_wes_no_reorder02_inbin(float2* d_kernel, float2* d_signal, float *d_ffdot_pw, int sigblock, int sig_tot_convlen, int sig_totlen, int offset, float scale, float2 *ip_edge_points)
/* convolution kernel using Karel Adamek's custom FFT, deployed here with modifications by Wes Armour.
   It performs the forward FFT and then loops through filters.(1-d blocks)
   It also uses half the templates and computes the rest using the complex conjugate.
   The modifications are optimizing speed by computing FFTs on two templates with one synchronization point. 
   Karel's kernel with removed de-shuffling of the data during the FFT.
   In this version we apply interbinning (fourier amplitudes interpolation in two bins)*/
{
  int tx = threadIdx.x;
  int bx = blockIdx.x;
  float2 ffdotcpx;
  float2 ffdotcpx_trans;
  float2 local_data; //fft'd data

  float one;
  float two;
  float three;
  float four;

  float2 A_DFT_value, B_DFT_value, sA_DFT_value, sB_DFT_value;
  float2 ftemp2, ftemp, stemp2, stemp;
  float2 W;
  int local_id, warp_id;
  int j, m_param;
  //int load_id, i, n;
  int parity, itemp;
  int A_read_index,B_read_index;
  //  int A_write_index,B_write_index;
  int PoT, PoTp1, q;
  
  // TODO: interpolation arrays size should be determined at runtime 
  // to save shared memory and conditionals
  //  extern 
  __shared__ float s_output_ip[2*KERNLEN]; 
  //extern 
  __shared__ float s_output_ip_trans[2*KERNLEN];
  //  s_output_ip_trans = &s_output_ip[2*sigblock];

  //complex input arrays
  __shared__ float2 s_input[KERNLEN]; //signal input data to FFT
  __shared__ float2 s_input_trans[KERNLEN]; //static allocation
 
  //setup ip arrays

  int index =  2*bx*sigblock + tx;
  int tidx =  bx * blockDim.x + tx;
  int sig_inbin_len = 2*sig_totlen;


  //load signal data to shared memory
  s_input[tx]=__ldg(&d_signal[tidx]);
  __syncthreads();


  //call custom device fft
  do_FFT_no_reorder(s_input,KERNLEN,NEXP);

  local_data.x = s_input[tx].x;
  local_data.y = s_input[tx].y;
  // __syncthreads();

  //complex multiplication power calculation loop over template columns
#pragma unroll 
  for (int i = 0; i < ZMAX/2; i++){
    //   __syncthreads();
    one = local_data.x*__ldg(&d_kernel[i*KERNLEN + tx].x);
    two = local_data.y*__ldg(&d_kernel[i*KERNLEN + tx].y);
    three = local_data.x*__ldg(&d_kernel[i*KERNLEN + tx].y);
    four = local_data.y*__ldg(&d_kernel[i*KERNLEN + tx].x);

    ffdotcpx.x = ((one - two) * scale);
    ffdotcpx.y = ((three + four) * scale);
    s_input[tx] = ffdotcpx;

    ffdotcpx_trans.x = ((one + two) * scale);
    ffdotcpx_trans.y = ((four -three) * scale);
    s_input_trans[tx] = ffdotcpx_trans;

    __syncthreads();

    //inverse fft
    local_id = threadIdx.x & (WARP - 1);
    warp_id = threadIdx.x/WARP;
    //-----> FFT
    //-->
    if(threadIdx.x<KERNLEN/2){
      PoT=1;
      PoTp1=2;	
      //--> First iteration
      itemp=local_id - (local_id&4294967294);
      parity=(1-itemp*2);
      
      A_DFT_value=s_input[local_id + warp_id*2*WARP];
      sA_DFT_value=s_input_trans[local_id + warp_id*2*WARP];
     
      B_DFT_value=s_input[local_id + warp_id*2*WARP + WARP];
      sB_DFT_value=s_input_trans[local_id + warp_id*2*WARP + WARP];

      //1st template
      A_DFT_value.x=parity*A_DFT_value.x+ __shfl(A_DFT_value.x,local_id+parity);
      A_DFT_value.y=parity*A_DFT_value.y+ __shfl(A_DFT_value.y,local_id+parity);
      B_DFT_value.x=parity*B_DFT_value.x+ __shfl(B_DFT_value.x,local_id+parity);
      B_DFT_value.y=parity*B_DFT_value.y+ __shfl(B_DFT_value.y,local_id+parity);
      //2nd template
      sA_DFT_value.x=parity*sA_DFT_value.x+ __shfl(sA_DFT_value.x,local_id+parity);
      sA_DFT_value.y=parity*sA_DFT_value.y+ __shfl(sA_DFT_value.y,local_id+parity);
      sB_DFT_value.x=parity*sB_DFT_value.x+ __shfl(sB_DFT_value.x,local_id+parity);
      sB_DFT_value.y=parity*sB_DFT_value.y+ __shfl(sB_DFT_value.y,local_id+parity);
      
      //--> First iteration end
		
      PoT=2;
      PoTp1=4;
      float fPoT=2.0f;
      for(q=2;q<6;q++){
	m_param = local_id & (PoTp1 - 1);
	itemp=__float2int_rz(__fdividef(__int2float_rz(m_param),fPoT));
	parity=(1-itemp*2);
	  
	W.x=cosf( 2.0f*3.141592654f*fdividef( (float) m_param, (float) (PoT*2)) );
	W.y=sinf( 2.0f*3.141592654f*fdividef( (float) m_param, (float) (PoT*2)) );
	  
	A_read_index=local_id+parity*itemp*PoT;
	B_read_index=local_id+(1-itemp)*PoT;
			
	//1st template
	ftemp2.x=__shfl(A_DFT_value.x,B_read_index);
	ftemp2.y=__shfl(A_DFT_value.y,B_read_index);				
	//2nd template
	stemp2.x=__shfl(sA_DFT_value.x,B_read_index);
	stemp2.y=__shfl(sA_DFT_value.y,B_read_index);
	
	//1st template
	A_DFT_value.x=__shfl(A_DFT_value.x,A_read_index) + W.x*ftemp2.x-W.y*ftemp2.y;
	A_DFT_value.y=__shfl(A_DFT_value.y,A_read_index) + W.x*ftemp2.y+W.y*ftemp2.x;
	//2nd template
	sA_DFT_value.x=__shfl(sA_DFT_value.x,A_read_index) + W.x*stemp2.x-W.y*stemp2.y;
	sA_DFT_value.y=__shfl(sA_DFT_value.y,A_read_index) + W.x*stemp2.y+W.y*stemp2.x;
	
	//1st template	  
	ftemp.x=__shfl(B_DFT_value.x,B_read_index);
	ftemp.y=__shfl(B_DFT_value.y,B_read_index);
	//2nd template	  
	stemp.x=__shfl(sB_DFT_value.x,B_read_index);
	stemp.y=__shfl(sB_DFT_value.y,B_read_index);

	//1st template
	B_DFT_value.x=__shfl(B_DFT_value.x,A_read_index) + W.x*ftemp.x - W.y*ftemp.y;
	B_DFT_value.y=__shfl(B_DFT_value.y,A_read_index) + W.x*ftemp.y + W.y*ftemp.x;		
	//2nd template
	sB_DFT_value.x=__shfl(sB_DFT_value.x,A_read_index) + W.x*stemp.x - W.y*stemp.y;
	sB_DFT_value.y=__shfl(sB_DFT_value.y,A_read_index) + W.x*stemp.y + W.y*stemp.x;		
	
	PoT=PoT<<1;
	fPoT=fPoT*2.0f;
	PoTp1=PoTp1<<1;
      }
      //1st template
      s_input[local_id + warp_id*2*WARP]=A_DFT_value;
      s_input[local_id + warp_id*2*WARP + WARP]=B_DFT_value;
      //2nd template
      s_input_trans[local_id + warp_id*2*WARP]=sA_DFT_value;
      s_input_trans[local_id + warp_id*2*WARP + WARP]=sB_DFT_value;
      
    }
    __syncthreads();
      
    for(q=6;q<=NEXP;q++){
      if(threadIdx.x<KERNLEN/2){
	m_param = threadIdx.x & (PoT - 1);
	j=threadIdx.x>>(q-1);
	  
	W=Get_W_value_inverse(PoTp1,m_param);

	A_read_index=j*PoTp1 + m_param;
	B_read_index=j*PoTp1 + m_param + PoT;
	  
	//1st template
	ftemp  = s_input[A_read_index];
	ftemp2 = s_input[B_read_index];
	//2nd template
	stemp  = s_input_trans[A_read_index];
	stemp2 = s_input_trans[B_read_index];
	
	//1st template
	A_DFT_value.x=ftemp.x + W.x*ftemp2.x - W.y*ftemp2.y;
	A_DFT_value.y=ftemp.y + W.x*ftemp2.y + W.y*ftemp2.x;
	//2nd template
	sA_DFT_value.x=stemp.x + W.x*stemp2.x - W.y*stemp2.y;
	sA_DFT_value.y=stemp.y + W.x*stemp2.y + W.y*stemp2.x;
	
	//1st template	  
	B_DFT_value.x=ftemp.x - W.x*ftemp2.x + W.y*ftemp2.y;
	B_DFT_value.y=ftemp.y - W.x*ftemp2.y - W.y*ftemp2.x;
	//2nd template	  
	sB_DFT_value.x=stemp.x - W.x*stemp2.x + W.y*stemp2.y;
	sB_DFT_value.y=stemp.y - W.x*stemp2.y - W.y*stemp2.x;
	
	PoT=PoT<<1;
	PoTp1=PoTp1<<1;		
      }
      //__syncthreads();
      if(threadIdx.x<KERNLEN/2){
	//1st template	  
	s_input[A_read_index]=A_DFT_value;
	s_input[B_read_index]=B_DFT_value;
	//2nd template	  
	s_input_trans[A_read_index]=sA_DFT_value;
	s_input_trans[B_read_index]=sB_DFT_value;	
      }
      __syncthreads();
    }

    //interbinning
    if (tx==offset ){
      //1. write edge points of each result to global array for interpolation
      ip_edge_points[bx] = s_input[tx];
      ip_edge_points[bx + gridDim.x] = s_input_trans[tx];  
      
      //2. read them in the shared memory
      s_input[tx + sigblock ] =__ldg(&ip_edge_points[bx+1]);

      if( bx < gridDim.x -1)
	s_input_trans[tx + sigblock] =__ldg(&ip_edge_points[bx + gridDim.x + 1]);
    }
    __syncthreads();
  
    if(tx < sigblock){
      //existing points
      s_output_ip[2*tx] = pwcalc(s_input[tx +offset]);
      s_output_ip_trans[2*tx] = pwcalc(s_input_trans[tx +offset]);
      
      // interpolated points
      s_output_ip[2*tx+1] = 0.616850275f * ((s_input[tx+offset].x - s_input[tx+1 + offset].x) * (s_input[tx +offset].x - s_input[tx+1 +offset].x) + (s_input[tx+offset].y - s_input[tx+1+offset].y) * (s_input[tx+offset].y - s_input[tx+offset+1].y)); 
      s_output_ip_trans[2*tx+1] = 0.616850275 * ((s_input_trans[tx+offset].x - s_input_trans[tx+1+offset].x) * (s_input_trans[tx+offset].x - s_input_trans[tx+1+offset].x) + (s_input_trans[tx+offset].y - s_input_trans[tx+1+offset].y) * (s_input_trans[tx+offset].y - s_input_trans[tx+1+offset].y)); 
    }
    __syncthreads();
    
    if(tx < sigblock){
      // 1st template
      d_ffdot_pw[(i*sig_inbin_len + index)] = s_output_ip[tx];
      d_ffdot_pw[(i*sig_inbin_len + index + sigblock)] = s_output_ip[tx + sigblock];
      
      // 2nd tempalte
      d_ffdot_pw[(unsigned int)((ZMAX - i)*sig_inbin_len + index)] = s_output_ip_trans[tx];
      d_ffdot_pw[(unsigned int)((ZMAX - i)*sig_inbin_len + index + sigblock)] = s_output_ip_trans[tx + sigblock];
    }
  }
  //-------END OF SYMMETRIC TEMPLATES

  // now do z=0
  one = local_data.x*__ldg(&d_kernel[(ZMAX/2)*KERNLEN + tx].x);
  two = local_data.y*__ldg(&d_kernel[(ZMAX/2)*KERNLEN + tx].y);
  three = local_data.x*__ldg(&d_kernel[(ZMAX/2)*KERNLEN + tx].y);
  four = local_data.y*__ldg(&d_kernel[(ZMAX/2)*KERNLEN + tx].x);
  
  ffdotcpx.x = ((one - two) * scale);
  ffdotcpx.y = ((three + four) * scale);
  s_input[tx] = ffdotcpx;
  __syncthreads();
  //call inverse fft
  do_IFFT_no_reorder(s_input, KERNLEN ,NEXP);
  
  //do interbinning
  if (tx==offset ){
    //1. write edge points of each result to global array for interpolation
    ip_edge_points[bx] = s_input[tx];
      
    //2. read them in the shared memory
    s_input[tx + sigblock ] =__ldg(&ip_edge_points[bx+1]);
  }
  __syncthreads();
  
  if(tx < sigblock){
    s_output_ip[2*tx] = pwcalc(s_input[tx +offset]);

    s_output_ip[2*tx+1] = 0.616850275f * ((s_input[tx+offset].x - s_input[tx+1 + offset].x) * (s_input[tx +offset].x - s_input[tx+1 +offset].x) + (s_input[tx+offset].y - s_input[tx+1+offset].y) * (s_input[tx+offset].y - s_input[tx+offset+1].y)); 
  }
  // write powers 
  if(tx < sigblock){
    d_ffdot_pw[((ZMAX/2)*sig_inbin_len + index)] = s_output_ip[tx];
    d_ffdot_pw[((ZMAX/2)*sig_inbin_len + index + sigblock)] = s_output_ip[tx + sigblock];
  }
  //--------END
}

#endif
