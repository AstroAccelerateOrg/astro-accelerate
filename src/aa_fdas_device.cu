/* Device functions for acceleration search  */
/* 06/11/2014 Sofia Dimoudi sofia.dimoudi@oerc.ox.ac.uk */

#include "aa_fdas_device.hpp"
#include "aa_device_cuda_deprecated_wrappers.cuh"

namespace astroaccelerate {

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

  /** \brief Forward no reorder FFTs. */
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

  static __device__ __inline__ void do_FFT_CT_DIF_2elem_no_reorder(float2 *s_input){
    float2 A_DFT_value, B_DFT_value;
    float2 W;
    float2 Aftemp, Bftemp;

    int local_id, warp_id;
    int j, m_param, parity;
    int A_read_index, B_read_index;
    int PoT, PoTm1, q;
	
    local_id = threadIdx.x & (WARP - 1);
    warp_id = threadIdx.x/WARP;
	
	
    //-----> FFT
    //-->
    PoTm1 = (KERNLEN>>1);
    PoT   = KERNLEN;

    for(q=(NEXP-1);q>4;q--){
      __syncthreads();
      m_param = threadIdx.x & (PoTm1 - 1);
      j=threadIdx.x>>q;
		
      W=Get_W_value(PoT, m_param);

      A_read_index=j*PoT + m_param;
      B_read_index=j*PoT + m_param + PoTm1;
		
      Aftemp = s_input[A_read_index];
      Bftemp = s_input[B_read_index];
		
      A_DFT_value.x = Aftemp.x + Bftemp.x;
      A_DFT_value.y = Aftemp.y + Bftemp.y;
		
      B_DFT_value.x = W.x*(Aftemp.x - Bftemp.x) - W.y*(Aftemp.y - Bftemp.y);
      B_DFT_value.y = W.x*(Aftemp.y - Bftemp.y) + W.y*(Aftemp.x - Bftemp.x);
		
      s_input[A_read_index]=A_DFT_value;
      s_input[B_read_index]=B_DFT_value;
		
      PoT=PoT>>1;
      PoTm1=PoTm1>>1;
    }

    __syncthreads();
    A_DFT_value=s_input[local_id + warp_id*2*WARP];
    B_DFT_value=s_input[local_id + warp_id*2*WARP + WARP];
	
    for(q=4;q>=0;q--){
      m_param = (local_id & (PoT - 1));
      j = m_param>>q;
      parity=(1-j*2);
      W = Get_W_value(PoT, j*(m_param-PoTm1));
		
      Aftemp.x = parity*A_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value.x, PoTm1);
      Aftemp.y = parity*A_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value.y, PoTm1);
      Bftemp.x = parity*B_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value.x, PoTm1);
      Bftemp.y = parity*B_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value.y, PoTm1);
		
      A_DFT_value.x = W.x*Aftemp.x - W.y*Aftemp.y; 
      A_DFT_value.y = W.x*Aftemp.y + W.y*Aftemp.x;
      B_DFT_value.x = W.x*Bftemp.x - W.y*Bftemp.y; 
      B_DFT_value.y = W.x*Bftemp.y + W.y*Bftemp.x;
		
      PoT=PoT>>1;
      PoTm1=PoTm1>>1;
    }
	
    s_input[local_id + warp_id*2*WARP] = A_DFT_value;
    s_input[local_id + warp_id*2*WARP + WARP] = B_DFT_value;
	
    __syncthreads();
  }

  static __device__ __inline__ void do_FFT_CT_DIF_4elem_no_reorder(float2 *s_input){
    float2 A_DFT_value, B_DFT_value, C_DFT_value, D_DFT_value;
    float2 W;
    float2 Aftemp, Bftemp, Cftemp, Dftemp;

    int local_id, warp_id;
    int j, m_param, parity;
    int A_read_index, B_read_index, C_read_index, D_read_index;
    int PoT, PoTm1, q;
	
    local_id = threadIdx.x & (WARP - 1);
    warp_id = threadIdx.x/WARP;
	
	
    //-----> FFT
    //-->
    PoTm1 = (KERNLEN>>1);
    PoT   = KERNLEN;
	
    //Highest iteration
    m_param = threadIdx.x;
    j=0;
    A_read_index = m_param;
    B_read_index = m_param + PoTm1;
    C_read_index = m_param + (PoTm1>>1);
    D_read_index = m_param + 3*(PoTm1>>1);
	
    W=Get_W_value(PoT, m_param);
	
    Aftemp = s_input[A_read_index];
    Bftemp = s_input[B_read_index];
    Cftemp = s_input[C_read_index];
    Dftemp = s_input[D_read_index];
	
    A_DFT_value.x = Aftemp.x + Bftemp.x;
    A_DFT_value.y = Aftemp.y + Bftemp.y;
    B_DFT_value.x = W.x*(Aftemp.x - Bftemp.x) - W.y*(Aftemp.y - Bftemp.y);
    B_DFT_value.y = W.x*(Aftemp.y - Bftemp.y) + W.y*(Aftemp.x - Bftemp.x);
	
    C_DFT_value.x = Cftemp.x + Dftemp.x;
    C_DFT_value.y = Cftemp.y + Dftemp.y;
    D_DFT_value.x = W.y*(Cftemp.x - Dftemp.x) + W.x*(Cftemp.y - Dftemp.y);
    D_DFT_value.y = W.y*(Cftemp.y - Dftemp.y) - W.x*(Cftemp.x - Dftemp.x);
	
    s_input[A_read_index]=A_DFT_value;
    s_input[B_read_index]=B_DFT_value;
    s_input[C_read_index]=C_DFT_value;
    s_input[D_read_index]=D_DFT_value;
	
    PoT=PoT>>1;
    PoTm1=PoTm1>>1;
	
    for(q=(NEXP-2);q>4;q--){
      __syncthreads();
      m_param = threadIdx.x & (PoTm1 - 1);
      j=threadIdx.x>>q;
		
      W=Get_W_value(PoT, m_param);

      A_read_index=j*(PoT<<1) + m_param;
      B_read_index=j*(PoT<<1) + m_param + PoTm1;
      C_read_index=j*(PoT<<1) + m_param + PoT;
      D_read_index=j*(PoT<<1) + m_param + 3*PoTm1;
		
      Aftemp = s_input[A_read_index];
      Bftemp = s_input[B_read_index];
      Cftemp = s_input[C_read_index];
      Dftemp = s_input[D_read_index];
		
      A_DFT_value.x = Aftemp.x + Bftemp.x;
      A_DFT_value.y = Aftemp.y + Bftemp.y;
      C_DFT_value.x = Cftemp.x + Dftemp.x;
      C_DFT_value.y = Cftemp.y + Dftemp.y;
		
      B_DFT_value.x = W.x*(Aftemp.x - Bftemp.x) - W.y*(Aftemp.y - Bftemp.y);
      B_DFT_value.y = W.x*(Aftemp.y - Bftemp.y) + W.y*(Aftemp.x - Bftemp.x);
      D_DFT_value.x = W.x*(Cftemp.x - Dftemp.x) - W.y*(Cftemp.y - Dftemp.y);
      D_DFT_value.y = W.x*(Cftemp.y - Dftemp.y) + W.y*(Cftemp.x - Dftemp.x);
		
      s_input[A_read_index]=A_DFT_value;
      s_input[B_read_index]=B_DFT_value;
      s_input[C_read_index]=C_DFT_value;
      s_input[D_read_index]=D_DFT_value;
		
      PoT=PoT>>1;
      PoTm1=PoTm1>>1;
    }

    __syncthreads();
    j = local_id + (warp_id<<2)*WARP;
    A_DFT_value = s_input[j];
    B_DFT_value = s_input[j + WARP];
    C_DFT_value = s_input[j + 2*WARP];
    D_DFT_value = s_input[j + 3*WARP];
	
    for(q=4;q>=0;q--){
      m_param = (local_id & (PoT - 1));
      j = m_param>>q;
      parity=(1-j*2);
      W = Get_W_value(PoT, j*(m_param-PoTm1));
		
      Aftemp.x = parity*A_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value.x, PoTm1);
      Aftemp.y = parity*A_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value.y, PoTm1);
      Bftemp.x = parity*B_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value.x, PoTm1);
      Bftemp.y = parity*B_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value.y, PoTm1);
      Cftemp.x = parity*C_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,C_DFT_value.x, PoTm1);
      Cftemp.y = parity*C_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,C_DFT_value.y, PoTm1);
      Dftemp.x = parity*D_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,D_DFT_value.x, PoTm1);
      Dftemp.y = parity*D_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,D_DFT_value.y, PoTm1);
		
      A_DFT_value.x = W.x*Aftemp.x - W.y*Aftemp.y; 
      A_DFT_value.y = W.x*Aftemp.y + W.y*Aftemp.x;
      B_DFT_value.x = W.x*Bftemp.x - W.y*Bftemp.y; 
      B_DFT_value.y = W.x*Bftemp.y + W.y*Bftemp.x;
      C_DFT_value.x = W.x*Cftemp.x - W.y*Cftemp.y; 
      C_DFT_value.y = W.x*Cftemp.y + W.y*Cftemp.x;
      D_DFT_value.x = W.x*Dftemp.x - W.y*Dftemp.y; 
      D_DFT_value.y = W.x*Dftemp.y + W.y*Dftemp.x;
		
      PoT=PoT>>1;
      PoTm1=PoTm1>>1;
    }
	
    j = local_id + (warp_id<<2)*WARP;
    s_input[j]          = A_DFT_value;
    s_input[j + WARP]   = B_DFT_value;
    s_input[j + 2*WARP] = C_DFT_value;
    s_input[j + 3*WARP] = D_DFT_value;
	
    __syncthreads();
  }

  //----------------- Forward no reorder FFTs

  //----------------- Inverse no reorder FFTs

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
		
      A_DFT_value.x=parity*A_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,A_DFT_value.x,local_id+parity);
      A_DFT_value.y=parity*A_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,A_DFT_value.y,local_id+parity);
		
      B_DFT_value.x=parity*B_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,B_DFT_value.x,local_id+parity);
      B_DFT_value.y=parity*B_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,B_DFT_value.y,local_id+parity);
		
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
			
	ftemp2.x=aa_shfl(AA_ASSUME_MASK,A_DFT_value.x,B_read_index);
	ftemp2.y=aa_shfl(AA_ASSUME_MASK,A_DFT_value.y,B_read_index);					
	A_DFT_value.x=aa_shfl(AA_ASSUME_MASK,A_DFT_value.x,A_read_index) + W.x*ftemp2.x-W.y*ftemp2.y;
	A_DFT_value.y=aa_shfl(AA_ASSUME_MASK,A_DFT_value.y,A_read_index) + W.x*ftemp2.y+W.y*ftemp2.x;
			
	ftemp.x=aa_shfl(AA_ASSUME_MASK,B_DFT_value.x,B_read_index);
	ftemp.y=aa_shfl(AA_ASSUME_MASK,B_DFT_value.y,B_read_index);					
	B_DFT_value.x=aa_shfl(AA_ASSUME_MASK,B_DFT_value.x,A_read_index) + W.x*ftemp.x - W.y*ftemp.y;
	B_DFT_value.y=aa_shfl(AA_ASSUME_MASK,B_DFT_value.y,A_read_index) + W.x*ftemp.y + W.y*ftemp.x;		
			
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

  static __device__ __inline__ void do_IFFT_mk11_2elem_2vertical_no_reorder(float2 *s_input_1, float2 *s_input_2){
    float2 A_DFT_value1, B_DFT_value1;
    float2 A_DFT_value2, B_DFT_value2;
    float2 W;
    float2 Aftemp1, Bftemp1;
    float2 Aftemp2, Bftemp2;

    int local_id, warp_id;
    int j, m_param;
    int parity, itemp;
    int A_read_index, B_read_index;
    int PoT, PoTp1, q;
	
    local_id = threadIdx.x & (WARP - 1);
    warp_id = threadIdx.x/WARP;

    //-----> FFT
    //-->
    PoT=1;
    PoTp1=2;	

    //--> First iteration
    itemp=local_id&1;
    parity=(1-itemp*2);
    A_DFT_value1=s_input_1[local_id + (warp_id<<1)*WARP];
    B_DFT_value1=s_input_1[local_id + (warp_id<<1)*WARP + WARP];
    A_DFT_value2=s_input_2[local_id + (warp_id<<1)*WARP];
    B_DFT_value2=s_input_2[local_id + (warp_id<<1)*WARP + WARP];
	
    __syncthreads();
	
    A_DFT_value1.x=parity*A_DFT_value1.x + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value1.x,1);
    A_DFT_value1.y=parity*A_DFT_value1.y + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value1.y,1);
    B_DFT_value1.x=parity*B_DFT_value1.x + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value1.x,1);
    B_DFT_value1.y=parity*B_DFT_value1.y + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value1.y,1);
	
    A_DFT_value2.x=parity*A_DFT_value2.x + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value2.x,1);
    A_DFT_value2.y=parity*A_DFT_value2.y + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value2.y,1);
    B_DFT_value2.x=parity*B_DFT_value2.x + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value2.x,1);
    B_DFT_value2.y=parity*B_DFT_value2.y + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value2.y,1);
	
	
    //--> Second through Fifth iteration (no synchronization)
    PoT=2;
    PoTp1=4;
    for(q=1;q<5;q++){
      m_param = (local_id & (PoTp1 - 1));
      itemp = m_param>>q;
      parity=((itemp<<1)-1);
      W = Get_W_value_inverse(PoTp1, itemp*m_param);
		
      Aftemp1.x = W.x*A_DFT_value1.x - W.y*A_DFT_value1.y;
      Aftemp1.y = W.x*A_DFT_value1.y + W.y*A_DFT_value1.x;
      Bftemp1.x = W.x*B_DFT_value1.x - W.y*B_DFT_value1.y;
      Bftemp1.y = W.x*B_DFT_value1.y + W.y*B_DFT_value1.x;
		
      Aftemp2.x = W.x*A_DFT_value2.x - W.y*A_DFT_value2.y;
      Aftemp2.y = W.x*A_DFT_value2.y + W.y*A_DFT_value2.x;
      Bftemp2.x = W.x*B_DFT_value2.x - W.y*B_DFT_value2.y;
      Bftemp2.y = W.x*B_DFT_value2.y + W.y*B_DFT_value2.x;
		
      A_DFT_value1.x = Aftemp1.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp1.x,PoT);
      A_DFT_value1.y = Aftemp1.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp1.y,PoT);
      B_DFT_value1.x = Bftemp1.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp1.x,PoT);
      B_DFT_value1.y = Bftemp1.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp1.y,PoT);
		
      A_DFT_value2.x = Aftemp2.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp2.x,PoT);
      A_DFT_value2.y = Aftemp2.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp2.y,PoT);
      B_DFT_value2.x = Bftemp2.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp2.x,PoT);
      B_DFT_value2.y = Bftemp2.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp2.y,PoT);	
		
      PoT=PoT<<1;
      PoTp1=PoTp1<<1;
    }
	
    itemp = local_id + (warp_id<<1)*WARP;
    s_input_1[itemp]          = A_DFT_value1;
    s_input_1[itemp + WARP]   = B_DFT_value1;
    s_input_2[itemp]          = A_DFT_value2;
    s_input_2[itemp + WARP]   = B_DFT_value2;
	
    for(q=5;q<NEXP;q++){
      __syncthreads();
      m_param = threadIdx.x & (PoT - 1);
      j=threadIdx.x>>q;
		
      W=Get_W_value_inverse(PoTp1,m_param);

      A_read_index=j*PoTp1 + m_param;
      B_read_index=j*PoTp1 + m_param + PoT;
		
      Aftemp1 = s_input_1[A_read_index];
      Bftemp1 = s_input_1[B_read_index];
      A_DFT_value1.x=Aftemp1.x + W.x*Bftemp1.x - W.y*Bftemp1.y;
      A_DFT_value1.y=Aftemp1.y + W.x*Bftemp1.y + W.y*Bftemp1.x;		
      B_DFT_value1.x=Aftemp1.x - W.x*Bftemp1.x + W.y*Bftemp1.y;
      B_DFT_value1.y=Aftemp1.y - W.x*Bftemp1.y - W.y*Bftemp1.x;
		
      Aftemp2 = s_input_2[A_read_index];
      Bftemp2 = s_input_2[B_read_index];
      A_DFT_value2.x=Aftemp2.x + W.x*Bftemp2.x - W.y*Bftemp2.y;
      A_DFT_value2.y=Aftemp2.y + W.x*Bftemp2.y + W.y*Bftemp2.x;		
      B_DFT_value2.x=Aftemp2.x - W.x*Bftemp2.x + W.y*Bftemp2.y;
      B_DFT_value2.y=Aftemp2.y - W.x*Bftemp2.y - W.y*Bftemp2.x;
		
      s_input_1[A_read_index]=A_DFT_value1;
      s_input_1[B_read_index]=B_DFT_value1;
      s_input_2[A_read_index]=A_DFT_value2;
      s_input_2[B_read_index]=B_DFT_value2;
		
      PoT=PoT<<1;
      PoTp1=PoTp1<<1;
    }
	
    __syncthreads();
  }

  static __device__ __inline__ void do_IFFT_mk11_2elem_no_reorder(float2 *s_input){
    float2 A_DFT_value, B_DFT_value;
    float2 W;
    float2 ftemp, ftemp2;

    int local_id, warp_id;
    int j, m_param;
    int parity, itemp;
    int A_read_index,B_read_index;
    int PoT, PoTp1, q;
	
    local_id = threadIdx.x & (WARP - 1);
    warp_id = threadIdx.x/WARP;
	
    //-----> FFT
    //-->
    PoT=1;
    PoTp1=2;	

    //--> First iteration
    itemp=local_id&1;
    parity=(1-itemp*2);
    A_DFT_value=s_input[local_id + warp_id*2*WARP];
    B_DFT_value=s_input[local_id + warp_id*2*WARP + WARP];
	
    __syncthreads();
	
    A_DFT_value.x=parity*A_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value.x,1);
    A_DFT_value.y=parity*A_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value.y,1);
	
    B_DFT_value.x=parity*B_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value.x,1);
    B_DFT_value.y=parity*B_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value.y,1);
	
	
    //--> Second through Fifth iteration (no synchronization)
    PoT=2;
    PoTp1=4;
    for(q=1;q<5;q++){
      m_param = (local_id & (PoTp1 - 1));
      itemp = m_param>>q;
      parity=(itemp*2-1);
      W = Get_W_value_inverse(PoTp1, itemp*m_param);
		
      ftemp2.x = W.x*A_DFT_value.x - W.y*A_DFT_value.y;
      ftemp2.y = W.x*A_DFT_value.y + W.y*A_DFT_value.x;
      ftemp.x = W.x*B_DFT_value.x - W.y*B_DFT_value.y;
      ftemp.y = W.x*B_DFT_value.y + W.y*B_DFT_value.x;
		
      A_DFT_value.x = ftemp2.x + parity*aa_shfl_xor(AA_ASSUME_MASK,ftemp2.x,PoT);
      A_DFT_value.y = ftemp2.y + parity*aa_shfl_xor(AA_ASSUME_MASK,ftemp2.y,PoT);
      B_DFT_value.x = ftemp.x + parity*aa_shfl_xor(AA_ASSUME_MASK,ftemp.x,PoT);
      B_DFT_value.y = ftemp.y + parity*aa_shfl_xor(AA_ASSUME_MASK,ftemp.y,PoT);
		
      PoT=PoT<<1;
      PoTp1=PoTp1<<1;
    }
	
    s_input[local_id + warp_id*2*WARP]=A_DFT_value;
    s_input[local_id + warp_id*2*WARP + WARP]=B_DFT_value;
	
    for(q=5;q<NEXP;q++){
      __syncthreads();
      m_param = threadIdx.x & (PoT - 1);
      j=threadIdx.x>>q;
		
      W=Get_W_value_inverse(PoTp1,m_param);

      A_read_index=j*PoTp1 + m_param;
      B_read_index=j*PoTp1 + m_param + PoT;
		
      ftemp  = s_input[A_read_index];
      ftemp2 = s_input[B_read_index];
		
      A_DFT_value.x=ftemp.x + W.x*ftemp2.x - W.y*ftemp2.y;
      A_DFT_value.y=ftemp.y + W.x*ftemp2.y + W.y*ftemp2.x;
		
      B_DFT_value.x=ftemp.x - W.x*ftemp2.x + W.y*ftemp2.y;
      B_DFT_value.y=ftemp.y - W.x*ftemp2.y - W.y*ftemp2.x;
		
      s_input[A_read_index]=A_DFT_value;
      s_input[B_read_index]=B_DFT_value;
		
      PoT=PoT<<1;
      PoTp1=PoTp1<<1;
    }
	
    __syncthreads();
  }

  static __device__ __inline__ void do_IFFT_mk11_4elem_no_reorder(float2 *s_input){
    float2 A_DFT_value, B_DFT_value, C_DFT_value, D_DFT_value;
    float2 W;
    float2 Aftemp, Bftemp, Cftemp, Dftemp;

    int local_id, warp_id;
    int j, m_param;
    int parity, itemp;
    int A_read_index, B_read_index, C_read_index, D_read_index;
    int PoT, PoTp1, q;
	
    local_id = threadIdx.x & (WARP - 1);
    warp_id = threadIdx.x/WARP;
	
    //-----> FFT
    //-->
    PoT=1;
    PoTp1=2;	

    //--> First iteration
    itemp=local_id&1;
    parity=(1-itemp*2);
    A_DFT_value=s_input[local_id + (warp_id<<2)*WARP];
    B_DFT_value=s_input[local_id + (warp_id<<2)*WARP + WARP];
    C_DFT_value=s_input[local_id + (warp_id<<2)*WARP + 2*WARP];
    D_DFT_value=s_input[local_id + (warp_id<<2)*WARP + 3*WARP];
	
    __syncthreads();
	
    A_DFT_value.x=parity*A_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value.x,1);
    A_DFT_value.y=parity*A_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value.y,1);
    B_DFT_value.x=parity*B_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value.x,1);
    B_DFT_value.y=parity*B_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value.y,1);
    C_DFT_value.x=parity*C_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,C_DFT_value.x,1);
    C_DFT_value.y=parity*C_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,C_DFT_value.y,1);
    D_DFT_value.x=parity*D_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK,D_DFT_value.x,1);
    D_DFT_value.y=parity*D_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK,D_DFT_value.y,1);
	
	
    //--> Second through Fifth iteration (no synchronization)
    PoT=2;
    PoTp1=4;
    for(q=1;q<5;q++){
      m_param = (local_id & (PoTp1 - 1));
      itemp = m_param>>q;
      parity=((itemp<<1)-1);
      W = Get_W_value_inverse(PoTp1, itemp*m_param);
		
      Aftemp.x = W.x*A_DFT_value.x - W.y*A_DFT_value.y;
      Aftemp.y = W.x*A_DFT_value.y + W.y*A_DFT_value.x;
      Bftemp.x = W.x*B_DFT_value.x - W.y*B_DFT_value.y;
      Bftemp.y = W.x*B_DFT_value.y + W.y*B_DFT_value.x;
      Cftemp.x = W.x*C_DFT_value.x - W.y*C_DFT_value.y;
      Cftemp.y = W.x*C_DFT_value.y + W.y*C_DFT_value.x;
      Dftemp.x = W.x*D_DFT_value.x - W.y*D_DFT_value.y;
      Dftemp.y = W.x*D_DFT_value.y + W.y*D_DFT_value.x;
		
      A_DFT_value.x = Aftemp.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp.x,PoT);
      A_DFT_value.y = Aftemp.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp.y,PoT);
      B_DFT_value.x = Bftemp.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp.x,PoT);
      B_DFT_value.y = Bftemp.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp.y,PoT);
      C_DFT_value.x = Cftemp.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Cftemp.x,PoT);
      C_DFT_value.y = Cftemp.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Cftemp.y,PoT);
      D_DFT_value.x = Dftemp.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Dftemp.x,PoT);
      D_DFT_value.y = Dftemp.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Dftemp.y,PoT);	
		
      PoT=PoT<<1;
      PoTp1=PoTp1<<1;
    }
	
    itemp = local_id + (warp_id<<2)*WARP;
    s_input[itemp]          = A_DFT_value;
    s_input[itemp + WARP]   = B_DFT_value;
    s_input[itemp + 2*WARP] = C_DFT_value;
    s_input[itemp + 3*WARP] = D_DFT_value;
	
    for(q=5;q<(NEXP-1);q++){
      __syncthreads();
      m_param = threadIdx.x & (PoT - 1);
      j=threadIdx.x>>q;
		
      W=Get_W_value_inverse(PoTp1,m_param);

      A_read_index=j*(PoTp1<<1) + m_param;
      B_read_index=j*(PoTp1<<1) + m_param + PoT;
      C_read_index=j*(PoTp1<<1) + m_param + PoTp1;
      D_read_index=j*(PoTp1<<1) + m_param + 3*PoT;
		
      Aftemp = s_input[A_read_index];
      Bftemp = s_input[B_read_index];
      A_DFT_value.x=Aftemp.x + W.x*Bftemp.x - W.y*Bftemp.y;
      A_DFT_value.y=Aftemp.y + W.x*Bftemp.y + W.y*Bftemp.x;		
      B_DFT_value.x=Aftemp.x - W.x*Bftemp.x + W.y*Bftemp.y;
      B_DFT_value.y=Aftemp.y - W.x*Bftemp.y - W.y*Bftemp.x;
		
      Cftemp = s_input[C_read_index];
      Dftemp = s_input[D_read_index];
      C_DFT_value.x=Cftemp.x + W.x*Dftemp.x - W.y*Dftemp.y;
      C_DFT_value.y=Cftemp.y + W.x*Dftemp.y + W.y*Dftemp.x;		
      D_DFT_value.x=Cftemp.x - W.x*Dftemp.x + W.y*Dftemp.y;
      D_DFT_value.y=Cftemp.y - W.x*Dftemp.y - W.y*Dftemp.x;
		
      s_input[A_read_index]=A_DFT_value;
      s_input[B_read_index]=B_DFT_value;
      s_input[C_read_index]=C_DFT_value;
      s_input[D_read_index]=D_DFT_value;
		
      PoT=PoT<<1;
      PoTp1=PoTp1<<1;
    }
	
    //last iteration
    __syncthreads();
    m_param = threadIdx.x;
	
    W=Get_W_value_inverse(PoTp1,m_param);
    
    A_read_index = m_param;
    B_read_index = m_param + PoT;
    C_read_index = m_param + (PoT>>1);
    D_read_index = m_param + 3*(PoT>>1);
	
    Aftemp = s_input[A_read_index];
    Bftemp = s_input[B_read_index];
    A_DFT_value.x=Aftemp.x + W.x*Bftemp.x - W.y*Bftemp.y;
    A_DFT_value.y=Aftemp.y + W.x*Bftemp.y + W.y*Bftemp.x;		
    B_DFT_value.x=Aftemp.x - W.x*Bftemp.x + W.y*Bftemp.y;
    B_DFT_value.y=Aftemp.y - W.x*Bftemp.y - W.y*Bftemp.x;
	
    Cftemp = s_input[C_read_index];
    Dftemp = s_input[D_read_index];
    C_DFT_value.x=Cftemp.x - W.y*Dftemp.x - W.x*Dftemp.y;
    C_DFT_value.y=Cftemp.y - W.y*Dftemp.y + W.x*Dftemp.x;		
    D_DFT_value.x=Cftemp.x + W.y*Dftemp.x + W.x*Dftemp.y;
    D_DFT_value.y=Cftemp.y + W.y*Dftemp.y - W.x*Dftemp.x;
	
    s_input[A_read_index]=A_DFT_value;
    s_input[B_read_index]=B_DFT_value;
    s_input[C_read_index]=C_DFT_value;
    s_input[D_read_index]=D_DFT_value;

    __syncthreads();	
  }

  static __device__ __inline__ void do_IFFT_mk11_4elem_2vertical_no_reorder(float2 *s_input1, float2 *s_input2){
    float2 A_DFT_value1, B_DFT_value1, C_DFT_value1, D_DFT_value1;
    float2 A_DFT_value2, B_DFT_value2, C_DFT_value2, D_DFT_value2;
    float2 W;
    float2 Aftemp1, Bftemp1, Cftemp1, Dftemp1;
    float2 Aftemp2, Bftemp2, Cftemp2, Dftemp2;

    int local_id, warp_id;
    int j, m_param;
    int parity, itemp;
    int A_read_index, B_read_index, C_read_index, D_read_index;
    int PoT, PoTp1, q;
	
    local_id = threadIdx.x & (WARP - 1);
    warp_id = threadIdx.x/WARP;
	
    //-----> FFT
    //-->
    PoT=1;
    PoTp1=2;	

    //--> First iteration
    itemp=local_id&1;
    parity=(1-itemp*2);
    A_DFT_value1=s_input1[local_id + (warp_id<<2)*WARP];
    B_DFT_value1=s_input1[local_id + (warp_id<<2)*WARP + WARP];
    C_DFT_value1=s_input1[local_id + (warp_id<<2)*WARP + 2*WARP];
    D_DFT_value1=s_input1[local_id + (warp_id<<2)*WARP + 3*WARP];
    A_DFT_value2=s_input2[local_id + (warp_id<<2)*WARP];
    B_DFT_value2=s_input2[local_id + (warp_id<<2)*WARP + WARP];
    C_DFT_value2=s_input2[local_id + (warp_id<<2)*WARP + 2*WARP];
    D_DFT_value2=s_input2[local_id + (warp_id<<2)*WARP + 3*WARP];
	
    __syncthreads();
	
    A_DFT_value1.x=parity*A_DFT_value1.x + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value1.x,1);
    A_DFT_value1.y=parity*A_DFT_value1.y + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value1.y,1);
    B_DFT_value1.x=parity*B_DFT_value1.x + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value1.x,1);
    B_DFT_value1.y=parity*B_DFT_value1.y + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value1.y,1);
    C_DFT_value1.x=parity*C_DFT_value1.x + aa_shfl_xor(AA_ASSUME_MASK,C_DFT_value1.x,1);
    C_DFT_value1.y=parity*C_DFT_value1.y + aa_shfl_xor(AA_ASSUME_MASK,C_DFT_value1.y,1);
    D_DFT_value1.x=parity*D_DFT_value1.x + aa_shfl_xor(AA_ASSUME_MASK,D_DFT_value1.x,1);
    D_DFT_value1.y=parity*D_DFT_value1.y + aa_shfl_xor(AA_ASSUME_MASK,D_DFT_value1.y,1);
	
    A_DFT_value2.x=parity*A_DFT_value2.x + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value2.x,1);
    A_DFT_value2.y=parity*A_DFT_value2.y + aa_shfl_xor(AA_ASSUME_MASK,A_DFT_value2.y,1);
    B_DFT_value2.x=parity*B_DFT_value2.x + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value2.x,1);
    B_DFT_value2.y=parity*B_DFT_value2.y + aa_shfl_xor(AA_ASSUME_MASK,B_DFT_value2.y,1);
    C_DFT_value2.x=parity*C_DFT_value2.x + aa_shfl_xor(AA_ASSUME_MASK,C_DFT_value2.x,1);
    C_DFT_value2.y=parity*C_DFT_value2.y + aa_shfl_xor(AA_ASSUME_MASK,C_DFT_value2.y,1);
    D_DFT_value2.x=parity*D_DFT_value2.x + aa_shfl_xor(AA_ASSUME_MASK,D_DFT_value2.x,1);
    D_DFT_value2.y=parity*D_DFT_value2.y + aa_shfl_xor(AA_ASSUME_MASK,D_DFT_value2.y,1);
	
	
    //--> Second through Fifth iteration (no synchronization)
    PoT=2;
    PoTp1=4;
    for(q=1;q<5;q++){
      m_param = (local_id & (PoTp1 - 1));
      itemp = m_param>>q;
      parity=((itemp<<1)-1);
      W = Get_W_value_inverse(PoTp1, itemp*m_param);
		
      Aftemp1.x = W.x*A_DFT_value1.x - W.y*A_DFT_value1.y;
      Aftemp1.y = W.x*A_DFT_value1.y + W.y*A_DFT_value1.x;
      Bftemp1.x = W.x*B_DFT_value1.x - W.y*B_DFT_value1.y;
      Bftemp1.y = W.x*B_DFT_value1.y + W.y*B_DFT_value1.x;
      Cftemp1.x = W.x*C_DFT_value1.x - W.y*C_DFT_value1.y;
      Cftemp1.y = W.x*C_DFT_value1.y + W.y*C_DFT_value1.x;
      Dftemp1.x = W.x*D_DFT_value1.x - W.y*D_DFT_value1.y;
      Dftemp1.y = W.x*D_DFT_value1.y + W.y*D_DFT_value1.x;
		
      Aftemp2.x = W.x*A_DFT_value2.x - W.y*A_DFT_value2.y;
      Aftemp2.y = W.x*A_DFT_value2.y + W.y*A_DFT_value2.x;
      Bftemp2.x = W.x*B_DFT_value2.x - W.y*B_DFT_value2.y;
      Bftemp2.y = W.x*B_DFT_value2.y + W.y*B_DFT_value2.x;
      Cftemp2.x = W.x*C_DFT_value2.x - W.y*C_DFT_value2.y;
      Cftemp2.y = W.x*C_DFT_value2.y + W.y*C_DFT_value2.x;
      Dftemp2.x = W.x*D_DFT_value2.x - W.y*D_DFT_value2.y;
      Dftemp2.y = W.x*D_DFT_value2.y + W.y*D_DFT_value2.x;
		
      A_DFT_value1.x = Aftemp1.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp1.x,PoT);
      A_DFT_value1.y = Aftemp1.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp1.y,PoT);
      B_DFT_value1.x = Bftemp1.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp1.x,PoT);
      B_DFT_value1.y = Bftemp1.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp1.y,PoT);
      C_DFT_value1.x = Cftemp1.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Cftemp1.x,PoT);
      C_DFT_value1.y = Cftemp1.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Cftemp1.y,PoT);
      D_DFT_value1.x = Dftemp1.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Dftemp1.x,PoT);
      D_DFT_value1.y = Dftemp1.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Dftemp1.y,PoT);
		
      A_DFT_value2.x = Aftemp2.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp2.x,PoT);
      A_DFT_value2.y = Aftemp2.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Aftemp2.y,PoT);
      B_DFT_value2.x = Bftemp2.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp2.x,PoT);
      B_DFT_value2.y = Bftemp2.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Bftemp2.y,PoT);
      C_DFT_value2.x = Cftemp2.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Cftemp2.x,PoT);
      C_DFT_value2.y = Cftemp2.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Cftemp2.y,PoT);
      D_DFT_value2.x = Dftemp2.x + parity*aa_shfl_xor(AA_ASSUME_MASK,Dftemp2.x,PoT);
      D_DFT_value2.y = Dftemp2.y + parity*aa_shfl_xor(AA_ASSUME_MASK,Dftemp2.y,PoT);	
		
      PoT=PoT<<1;
      PoTp1=PoTp1<<1;
    }
	
    itemp = local_id + (warp_id<<2)*WARP;
    s_input1[itemp]          = A_DFT_value1;
    s_input1[itemp + WARP]   = B_DFT_value1;
    s_input1[itemp + 2*WARP] = C_DFT_value1;
    s_input1[itemp + 3*WARP] = D_DFT_value1;
	
    s_input2[itemp]          = A_DFT_value2;
    s_input2[itemp + WARP]   = B_DFT_value2;
    s_input2[itemp + 2*WARP] = C_DFT_value2;
    s_input2[itemp + 3*WARP] = D_DFT_value2;
	
    for(q=5;q<(NEXP-1);q++){
      __syncthreads();
      m_param = threadIdx.x & (PoT - 1);
      j=threadIdx.x>>q;
		
      W=Get_W_value_inverse(PoTp1,m_param);

      A_read_index=j*(PoTp1<<1) + m_param;
      B_read_index=j*(PoTp1<<1) + m_param + PoT;
      C_read_index=j*(PoTp1<<1) + m_param + PoTp1;
      D_read_index=j*(PoTp1<<1) + m_param + 3*PoT;
		
      Aftemp1 = s_input1[A_read_index];
      Bftemp1 = s_input1[B_read_index];
      A_DFT_value1.x=Aftemp1.x + W.x*Bftemp1.x - W.y*Bftemp1.y;
      A_DFT_value1.y=Aftemp1.y + W.x*Bftemp1.y + W.y*Bftemp1.x;		
      B_DFT_value1.x=Aftemp1.x - W.x*Bftemp1.x + W.y*Bftemp1.y;
      B_DFT_value1.y=Aftemp1.y - W.x*Bftemp1.y - W.y*Bftemp1.x;
		
      Aftemp2 = s_input2[A_read_index];
      Bftemp2 = s_input2[B_read_index];
      A_DFT_value2.x=Aftemp2.x + W.x*Bftemp2.x - W.y*Bftemp2.y;
      A_DFT_value2.y=Aftemp2.y + W.x*Bftemp2.y + W.y*Bftemp2.x;		
      B_DFT_value2.x=Aftemp2.x - W.x*Bftemp2.x + W.y*Bftemp2.y;
      B_DFT_value2.y=Aftemp2.y - W.x*Bftemp2.y - W.y*Bftemp2.x;
		
      Cftemp1 = s_input1[C_read_index];
      Dftemp1 = s_input1[D_read_index];
      C_DFT_value1.x=Cftemp1.x + W.x*Dftemp1.x - W.y*Dftemp1.y;
      C_DFT_value1.y=Cftemp1.y + W.x*Dftemp1.y + W.y*Dftemp1.x;		
      D_DFT_value1.x=Cftemp1.x - W.x*Dftemp1.x + W.y*Dftemp1.y;
      D_DFT_value1.y=Cftemp1.y - W.x*Dftemp1.y - W.y*Dftemp1.x;

      Cftemp2 = s_input2[C_read_index];
      Dftemp2 = s_input2[D_read_index];
      C_DFT_value2.x=Cftemp2.x + W.x*Dftemp2.x - W.y*Dftemp2.y;
      C_DFT_value2.y=Cftemp2.y + W.x*Dftemp2.y + W.y*Dftemp2.x;		
      D_DFT_value2.x=Cftemp2.x - W.x*Dftemp2.x + W.y*Dftemp2.y;
      D_DFT_value2.y=Cftemp2.y - W.x*Dftemp2.y - W.y*Dftemp2.x;
		
      s_input1[A_read_index]=A_DFT_value1;
      s_input1[B_read_index]=B_DFT_value1;
      s_input1[C_read_index]=C_DFT_value1;
      s_input1[D_read_index]=D_DFT_value1;
		
      s_input2[A_read_index]=A_DFT_value2;
      s_input2[B_read_index]=B_DFT_value2;
      s_input2[C_read_index]=C_DFT_value2;
      s_input2[D_read_index]=D_DFT_value2;
		
      PoT=PoT<<1;
      PoTp1=PoTp1<<1;
    }
	
    //last iteration
    __syncthreads();
    m_param = threadIdx.x;
	
    W=Get_W_value_inverse(PoTp1,m_param);
    
    A_read_index = m_param;
    B_read_index = m_param + PoT;
    C_read_index = m_param + (PoT>>1);
    D_read_index = m_param + 3*(PoT>>1);
	
    Aftemp1 = s_input1[A_read_index];
    Bftemp1 = s_input1[B_read_index];
    A_DFT_value1.x=Aftemp1.x + W.x*Bftemp1.x - W.y*Bftemp1.y;
    A_DFT_value1.y=Aftemp1.y + W.x*Bftemp1.y + W.y*Bftemp1.x;		
    B_DFT_value1.x=Aftemp1.x - W.x*Bftemp1.x + W.y*Bftemp1.y;
    B_DFT_value1.y=Aftemp1.y - W.x*Bftemp1.y - W.y*Bftemp1.x;

    Aftemp2 = s_input2[A_read_index];
    Bftemp2 = s_input2[B_read_index];
    A_DFT_value2.x=Aftemp2.x + W.x*Bftemp2.x - W.y*Bftemp2.y;
    A_DFT_value2.y=Aftemp2.y + W.x*Bftemp2.y + W.y*Bftemp2.x;		
    B_DFT_value2.x=Aftemp2.x - W.x*Bftemp2.x + W.y*Bftemp2.y;
    B_DFT_value2.y=Aftemp2.y - W.x*Bftemp2.y - W.y*Bftemp2.x;	
	
    Cftemp1 = s_input1[C_read_index];
    Dftemp1 = s_input1[D_read_index];
    C_DFT_value1.x=Cftemp1.x - W.y*Dftemp1.x - W.x*Dftemp1.y;
    C_DFT_value1.y=Cftemp1.y - W.y*Dftemp1.y + W.x*Dftemp1.x;		
    D_DFT_value1.x=Cftemp1.x + W.y*Dftemp1.x + W.x*Dftemp1.y;
    D_DFT_value1.y=Cftemp1.y + W.y*Dftemp1.y - W.x*Dftemp1.x;
	
    Cftemp2 = s_input2[C_read_index];
    Dftemp2 = s_input2[D_read_index];
    C_DFT_value2.x=Cftemp2.x - W.y*Dftemp2.x - W.x*Dftemp2.y;
    C_DFT_value2.y=Cftemp2.y - W.y*Dftemp2.y + W.x*Dftemp2.x;		
    D_DFT_value2.x=Cftemp2.x + W.y*Dftemp2.x + W.x*Dftemp2.y;
    D_DFT_value2.y=Cftemp2.y + W.y*Dftemp2.y - W.x*Dftemp2.x;
	
    s_input1[A_read_index]=A_DFT_value1;
    s_input1[B_read_index]=B_DFT_value1;
    s_input1[C_read_index]=C_DFT_value1;
    s_input1[D_read_index]=D_DFT_value1;	
	
    s_input2[A_read_index]=A_DFT_value2;
    s_input2[B_read_index]=B_DFT_value2;
    s_input2[C_read_index]=C_DFT_value2;
    s_input2[D_read_index]=D_DFT_value2;

    __syncthreads();	
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

  void call_kernel_cuda_overlap_copy(float2 *const d_ext_data, float2 *const d_cpx_signal, const int &sigblock, const int &sig_rfftlen, const int &sig_tot_convlen, const int &kern_offset, const int &total_blocks) {
    cuda_overlap_copy<<<KERNLEN/64, 64 >>>(d_ext_data, d_cpx_signal, sigblock, sig_rfftlen, sig_tot_convlen, kern_offset, total_blocks);
  }

  __global__ void cuda_overlap_copy_smallblk(float2* d_ext_data, float2* d_cpx_signal, int sigblock, int sig_rfftlen, int sig_tot_convlen, int kern_offset, int total_blocks)
  {
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int read_idx = blockIdx.x*sigblock - kern_offset + threadIdx.x;
    int write_idx = blockIdx.x*KERNLEN + threadIdx.x;

    //initialize the array  
    if (tid < sig_tot_convlen){
      d_ext_data[tid].x = 0.0f;
      d_ext_data[tid].y = 0.0f;
    }

    if (threadIdx.x >= kern_offset && blockIdx.x == 0 ) //first block
      d_ext_data[threadIdx.x] = d_cpx_signal[threadIdx.x - kern_offset ];
  
    // copy overlapped blocks
    if (blockIdx.x > 0 && read_idx < sig_rfftlen){    
      d_ext_data[write_idx] = d_cpx_signal[read_idx];
    }
  }

  void call_kernel_cuda_overlap_copy_smallblk(const int &blocks, float2 *const d_ext_data, float2 *const d_cpx_signal, const int &sigblock, const int &sig_rfftlen, const int &sig_tot_convlen, const int &kern_offset, const int &total_blocks) {
    cuda_overlap_copy_smallblk<<<blocks, KERNLEN>>>(d_ext_data, d_cpx_signal, sigblock, sig_rfftlen, sig_tot_convlen, kern_offset, total_blocks);
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

  void call_kernel_cuda_convolve_reg_1d_halftemps(const int &blocks, const int &threads, float2 *const d_kernel, float2 *const d_signal, float2 *const d_ffdot_plane, const int &sig_tot_convlen, const float &scale) {
    cuda_convolve_reg_1d_halftemps<<<blocks, threads>>>(d_kernel, d_signal, d_ffdot_plane, sig_tot_convlen, scale);
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

  void call_kernel_cuda_ffdotpow_concat_2d(const dim3 &blocks, const dim3 &threads, float2 *const d_ffdot_plane_cpx, float *const d_ffdot_plane, const int &sigblock, const int &kern_offset, const int &total_blocks, const int &sig_tot_convlen, const int &sig_totlen) {
    cuda_ffdotpow_concat_2d<<<blocks, threads>>>(d_ffdot_plane_cpx, d_ffdot_plane, sigblock, kern_offset, total_blocks, sig_tot_convlen, sig_totlen);
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

  void call_kernel_cuda_ffdotpow_concat_2d_inbin(const dim3 &blocks, const dim3 &threads, float2 *const d_ffdot_plane_cpx, float *const d_ffdot_plane, const int &sigblock, const int &kern_offset, const int &total_blocks, const int &sig_tot_convlen, const int &sig_totlen) {
    cuda_ffdotpow_concat_2d_inbin<<<blocks, threads>>>(d_ffdot_plane_cpx, d_ffdot_plane, sigblock, kern_offset, total_blocks, sig_tot_convlen, sig_totlen);
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

  void call_kernel_customfft_fwd_temps_no_reorder(float2 *const d_signal) {
    customfft_fwd_temps_no_reorder<<<NKERN,KERNLEN>>>( d_signal);
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
	A_DFT_value.x=parity*A_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,A_DFT_value.x,local_id+parity);
	A_DFT_value.y=parity*A_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,A_DFT_value.y,local_id+parity);
	B_DFT_value.x=parity*B_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,B_DFT_value.x,local_id+parity);
	B_DFT_value.y=parity*B_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,B_DFT_value.y,local_id+parity);
	//2nd template

	sA_DFT_value.x=parity*sA_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,sA_DFT_value.x,local_id+parity);
	sA_DFT_value.y=parity*sA_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,sA_DFT_value.y,local_id+parity);
	sB_DFT_value.x=parity*sB_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,sB_DFT_value.x,local_id+parity);
	sB_DFT_value.y=parity*sB_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,sB_DFT_value.y,local_id+parity);

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
	  ftemp2.x=aa_shfl(AA_ASSUME_MASK,A_DFT_value.x,B_read_index);
	  ftemp2.y=aa_shfl(AA_ASSUME_MASK,A_DFT_value.y,B_read_index);
	  //2nd template
	  stemp2.x=aa_shfl(AA_ASSUME_MASK,sA_DFT_value.x,B_read_index);
	  stemp2.y=aa_shfl(AA_ASSUME_MASK,sA_DFT_value.y,B_read_index);
	
	  //1st template
	  A_DFT_value.x=aa_shfl(AA_ASSUME_MASK,A_DFT_value.x,A_read_index) + W.x*ftemp2.x-W.y*ftemp2.y;
	  A_DFT_value.y=aa_shfl(AA_ASSUME_MASK,A_DFT_value.y,A_read_index) + W.x*ftemp2.y+W.y*ftemp2.x;
	  //2nd template
	  sA_DFT_value.x=aa_shfl(AA_ASSUME_MASK,sA_DFT_value.x,A_read_index) + W.x*stemp2.x-W.y*stemp2.y;
	  sA_DFT_value.y=aa_shfl(AA_ASSUME_MASK,sA_DFT_value.y,A_read_index) + W.x*stemp2.y+W.y*stemp2.x;
	

	  //1st template	  
	  ftemp.x=aa_shfl(AA_ASSUME_MASK,B_DFT_value.x,B_read_index);
	  ftemp.y=aa_shfl(AA_ASSUME_MASK,B_DFT_value.y,B_read_index);
	  //2nd template	  
	  stemp.x=aa_shfl(AA_ASSUME_MASK,sB_DFT_value.x,B_read_index);
	  stemp.y=aa_shfl(AA_ASSUME_MASK,sB_DFT_value.y,B_read_index);
	
	  //1st template
	  B_DFT_value.x=aa_shfl(AA_ASSUME_MASK,B_DFT_value.x,A_read_index) + W.x*ftemp.x - W.y*ftemp.y;
	  B_DFT_value.y=aa_shfl(AA_ASSUME_MASK,B_DFT_value.y,A_read_index) + W.x*ftemp.y + W.y*ftemp.x;
	  //2nd template
	  sB_DFT_value.x=aa_shfl(AA_ASSUME_MASK,sB_DFT_value.x,A_read_index) + W.x*stemp.x - W.y*stemp.y;
	  sB_DFT_value.y=aa_shfl(AA_ASSUME_MASK,sB_DFT_value.y,A_read_index) + W.x*stemp.y + W.y*stemp.x;
	  
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
	__syncthreads();
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

	  //1st template	  
	  s_input[A_read_index]=A_DFT_value;
	  s_input[B_read_index]=B_DFT_value;
	  //2nd template	  
	  s_input_trans[A_read_index]=sA_DFT_value;
	  s_input_trans[B_read_index]=sB_DFT_value;
	}
      }
        
      __syncthreads();
      if(tx >= offset && tx < sigblock + offset){
	d_ffdot_pw[i*sig_totlen + index - offset] = pwcalc(s_input[tx]);

	d_ffdot_pw[(ZMAX - i )*sig_totlen + index - offset] = pwcalc(s_input_trans[tx]);
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

  void call_kernel_cuda_convolve_customfft_wes_no_reorder02(const int &blocks, float2 *const d_kernel, float2 *const d_signal, float *const d_ffdot_pw, const int &sigblock, const int &sig_tot_convlen, const int &sig_totlen, const int &offset, const float &scale) {
    cuda_convolve_customfft_wes_no_reorder02<<<blocks, KERNLEN>>>(d_kernel, d_signal, d_ffdot_pw, sigblock, sig_tot_convlen, sig_totlen, offset, scale); 
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
	A_DFT_value.x=parity*A_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,A_DFT_value.x,local_id+parity);
	A_DFT_value.y=parity*A_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,A_DFT_value.y,local_id+parity);
	B_DFT_value.x=parity*B_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,B_DFT_value.x,local_id+parity);
	B_DFT_value.y=parity*B_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,B_DFT_value.y,local_id+parity);
	//2nd template
	sA_DFT_value.x=parity*sA_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,sA_DFT_value.x,local_id+parity);
	sA_DFT_value.y=parity*sA_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,sA_DFT_value.y,local_id+parity);
	sB_DFT_value.x=parity*sB_DFT_value.x+ aa_shfl(AA_ASSUME_MASK,sB_DFT_value.x,local_id+parity);
	sB_DFT_value.y=parity*sB_DFT_value.y+ aa_shfl(AA_ASSUME_MASK,sB_DFT_value.y,local_id+parity);
      
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
	  ftemp2.x=aa_shfl(AA_ASSUME_MASK,A_DFT_value.x,B_read_index);
	  ftemp2.y=aa_shfl(AA_ASSUME_MASK,A_DFT_value.y,B_read_index);
	  //2nd template
	  stemp2.x=aa_shfl(AA_ASSUME_MASK,sA_DFT_value.x,B_read_index);
	  stemp2.y=aa_shfl(AA_ASSUME_MASK,sA_DFT_value.y,B_read_index);
	
	  //1st template
	  A_DFT_value.x=aa_shfl(AA_ASSUME_MASK,A_DFT_value.x,A_read_index) + W.x*ftemp2.x-W.y*ftemp2.y;
	  A_DFT_value.y=aa_shfl(AA_ASSUME_MASK,A_DFT_value.y,A_read_index) + W.x*ftemp2.y+W.y*ftemp2.x;
	  //2nd template
	  sA_DFT_value.x=aa_shfl(AA_ASSUME_MASK,sA_DFT_value.x,A_read_index) + W.x*stemp2.x-W.y*stemp2.y;
	  sA_DFT_value.y=aa_shfl(AA_ASSUME_MASK,sA_DFT_value.y,A_read_index) + W.x*stemp2.y+W.y*stemp2.x;
	
	  //1st template	  
	  ftemp.x=aa_shfl(AA_ASSUME_MASK,B_DFT_value.x,B_read_index);
	  ftemp.y=aa_shfl(AA_ASSUME_MASK,B_DFT_value.y,B_read_index);
	  //2nd template	  
	  stemp.x=aa_shfl(AA_ASSUME_MASK,sB_DFT_value.x,B_read_index);
	  stemp.y=aa_shfl(AA_ASSUME_MASK,sB_DFT_value.y,B_read_index);

	  //1st template
	  B_DFT_value.x=aa_shfl(AA_ASSUME_MASK,B_DFT_value.x,A_read_index) + W.x*ftemp.x - W.y*ftemp.y;
	  B_DFT_value.y=aa_shfl(AA_ASSUME_MASK,B_DFT_value.y,A_read_index) + W.x*ftemp.y + W.y*ftemp.x;
	  //2nd template
	  sB_DFT_value.x=aa_shfl(AA_ASSUME_MASK,sB_DFT_value.x,A_read_index) + W.x*stemp.x - W.y*stemp.y;
	  sB_DFT_value.y=aa_shfl(AA_ASSUME_MASK,sB_DFT_value.y,A_read_index) + W.x*stemp.y + W.y*stemp.x;
	
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

  void call_kernel_cuda_convolve_customfft_wes_no_reorder02_inbin(const int &blocks, float2 *const d_kernel, float2 *const d_signal, float *const d_ffdot_pw, const int &sigblock, const int &sig_tot_convlen, const int &sig_totlen, const int &offset, const float &scale, float2 *const ip_edge_points) {
    cuda_convolve_customfft_wes_no_reorder02_inbin<<<blocks, KERNLEN>>>(d_kernel, d_signal, d_ffdot_pw, sigblock, sig_tot_convlen, sig_totlen, offset, scale, ip_edge_points);
  }


  __global__ void GPU_CONV_kFFT_mk11_2elem_2v(float2 const* __restrict__ d_input_signal, float *d_output_plane_reduced, float2 const* __restrict__ d_templates, int useful_part_size, int offset, int nConvolutions, float scale) {
    __shared__ float2 s_input_1[KERNLEN];
    __shared__ float2 s_input_2[KERNLEN];
    // Convolution
    float2 r_templates_1[2];
    float2 r_templates_2[2];
    float2 signal[2];
    int pos, t;
    // Loading data
    //prepare_signal(s_input_1, d_input_signal, useful_part_size, offset);
    s_input_1[threadIdx.x]=__ldg(&d_input_signal[blockIdx.x*KERNLEN + threadIdx.x]);
    s_input_1[threadIdx.x + (KERNLEN>>1)]=__ldg(&d_input_signal[blockIdx.x*KERNLEN + threadIdx.x + (KERNLEN>>1)]);

    do_FFT_CT_DIF_2elem_no_reorder(s_input_1);
	
    signal[0]=s_input_1[threadIdx.x];
    signal[1]=s_input_1[threadIdx.x + (KERNLEN>>1)];
	
    for(t=0; t<ZMAX/2; t++){
      // Loading templates
      pos = t*KERNLEN + threadIdx.x;
      r_templates_1[0]=__ldg(&d_templates[pos]);
      r_templates_1[1]=__ldg(&d_templates[pos + (KERNLEN>>1)]);
      r_templates_2[0]=__ldg(&d_templates[pos]);
      r_templates_2[1]=__ldg(&d_templates[pos + (KERNLEN>>1)]);

      // Convolution
      s_input_1[threadIdx.x].x = (r_templates_1[0].x*signal[0].x - r_templates_1[0].y*signal[0].y)*scale;
      s_input_1[threadIdx.x].y = (r_templates_1[0].x*signal[0].y + r_templates_1[0].y*signal[0].x)*scale;
		
      s_input_2[threadIdx.x].x = (r_templates_2[0].x*signal[0].x + r_templates_2[0].y*signal[0].y)*scale;
      s_input_2[threadIdx.x].y = (r_templates_2[0].x*signal[0].y - r_templates_2[0].y*signal[0].x)*scale;
		
      s_input_1[threadIdx.x + (KERNLEN>>1)].x = (r_templates_1[1].x*signal[1].x - r_templates_1[1].y*signal[1].y)*scale;
      s_input_1[threadIdx.x + (KERNLEN>>1)].y = (r_templates_1[1].x*signal[1].y + r_templates_1[1].y*signal[1].x)*scale;
		
      s_input_2[threadIdx.x + (KERNLEN>>1)].x = (r_templates_2[1].x*signal[1].x + r_templates_2[1].y*signal[1].y)*scale;
      s_input_2[threadIdx.x + (KERNLEN>>1)].y = (r_templates_2[1].x*signal[1].y - r_templates_2[1].y*signal[1].x)*scale;
		
      __syncthreads();
		
      //----------> IFFT
      do_IFFT_mk11_2elem_2vertical_no_reorder(s_input_1, s_input_2);
      //----------<
		
      // Saving data
      pos = t*useful_part_size*nConvolutions + blockIdx.x*useful_part_size + threadIdx.x;
      if( threadIdx.x>=offset && threadIdx.x<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos - offset] = pwcalc(s_input_1[threadIdx.x]);
      }
      if( (threadIdx.x+(KERNLEN>>1))>=offset && (threadIdx.x+(KERNLEN>>1))<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos + (KERNLEN>>1) - offset] = pwcalc(s_input_1[threadIdx.x + (KERNLEN>>1)]);
      }
      pos = (ZMAX-t)*useful_part_size*nConvolutions + blockIdx.x*useful_part_size + threadIdx.x;
      if( threadIdx.x>=offset && threadIdx.x<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos - offset] = pwcalc(s_input_2[threadIdx.x]);
      }
      if( (threadIdx.x+(KERNLEN>>1))>=offset && (threadIdx.x+(KERNLEN>>1))<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos + (KERNLEN>>1) - offset] = pwcalc(s_input_2[threadIdx.x + (KERNLEN>>1)]);
      }
    }
	
    // now do z=0
    pos = (ZMAX>>1)*KERNLEN + threadIdx.x;
    r_templates_1[0]=__ldg(&d_templates[pos]);
    r_templates_1[1]=__ldg(&d_templates[pos + (KERNLEN>>1)]);
	
    s_input_1[threadIdx.x].x=(r_templates_1[0].x*signal[0].x - r_templates_1[0].y*signal[0].y)*scale;
    s_input_1[threadIdx.x].y=(r_templates_1[0].x*signal[0].y + r_templates_1[0].y*signal[0].x)*scale;
    s_input_1[threadIdx.x + (KERNLEN>>1)].x=(r_templates_1[1].x*signal[1].x - r_templates_1[1].y*signal[1].y)*scale;
    s_input_1[threadIdx.x + (KERNLEN>>1)].y=(r_templates_1[1].x*signal[1].y + r_templates_1[1].y*signal[1].x)*scale;
	
    __syncthreads();
	
    //call inverse fft
    do_IFFT_mk11_2elem_no_reorder(s_input_1);
	
    pos = (ZMAX/2)*useful_part_size*nConvolutions + blockIdx.x*useful_part_size + threadIdx.x;
    if( threadIdx.x>=offset && threadIdx.x<(useful_part_size+offset) ) {
      d_output_plane_reduced[pos - offset] = pwcalc(s_input_1[threadIdx.x]);
    }
    if( (threadIdx.x+(KERNLEN>>1))>=offset && (threadIdx.x+(KERNLEN>>1))<(useful_part_size+offset) ) {
      d_output_plane_reduced[pos + (KERNLEN>>1) - offset] = pwcalc(s_input_1[threadIdx.x + (KERNLEN>>1)]);
    }
  }


  __global__ void GPU_CONV_kFFT_mk11_4elem_2v(float2 const* __restrict__ d_input_signal, float *d_output_plane_reduced, float2 const* __restrict__ d_templates, int useful_part_size, int offset, int nConvolutions, float scale) {
    __shared__ float2 s_input_1[KERNLEN];
    __shared__ float2 s_input_2[KERNLEN];
    // Convolution
    float2 r_templates_1[4];
    //float2 r_templates_2[4];
    float2 signal[4];
    float xx, yy, xy, yx;
    int pos, t;
    // Loading data
    //prepare_signal(s_input_1, d_input_signal, useful_part_size, offset);
	
    s_input_1[threadIdx.x]                  = __ldg(&d_input_signal[blockIdx.x*KERNLEN + threadIdx.x]);
    s_input_1[threadIdx.x + (KERNLEN>>2)]   = __ldg(&d_input_signal[blockIdx.x*KERNLEN + threadIdx.x + (KERNLEN>>2)]);
    s_input_1[threadIdx.x + (KERNLEN>>1)]   = __ldg(&d_input_signal[blockIdx.x*KERNLEN + threadIdx.x + (KERNLEN>>1)]);
    s_input_1[threadIdx.x + 3*(KERNLEN>>2)] = __ldg(&d_input_signal[blockIdx.x*KERNLEN + threadIdx.x + 3*(KERNLEN>>2)]);

    do_FFT_CT_DIF_4elem_no_reorder(s_input_1);
	
    signal[0] = s_input_1[threadIdx.x];
    signal[1] = s_input_1[threadIdx.x + (KERNLEN>>2)];
    signal[2] = s_input_1[threadIdx.x + (KERNLEN>>1)];
    signal[3] = s_input_1[threadIdx.x + 3*(KERNLEN>>2)];
	
    d_output_plane_reduced[threadIdx.x] = s_input_1[threadIdx.x].x;
	
	
    for(t=0; t<ZMAX/2; t++){
      // Loading templates
      pos = t*KERNLEN + threadIdx.x;
      r_templates_1[0] = __ldg(&d_templates[pos]);
      r_templates_1[1] = __ldg(&d_templates[pos + (KERNLEN>>2)]);
      r_templates_1[2] = __ldg(&d_templates[pos + (KERNLEN>>1)]);
      r_templates_1[3] = __ldg(&d_templates[pos + 3*(KERNLEN>>2)]);
		
      //r_templates_2[0] = __ldg(&d_templates[pos]);
      //r_templates_2[1] = __ldg(&d_templates[pos + (KERNLEN>>2)]);
      //r_templates_2[2] = __ldg(&d_templates[pos + (KERNLEN>>1)]);
      //r_templates_2[3] = __ldg(&d_templates[pos + 3*(KERNLEN>>2)]);
		
	
      // Convolution
      xx = r_templates_1[0].x*signal[0].x;
      yy = r_templates_1[0].y*signal[0].y;
      xy = r_templates_1[0].x*signal[0].y;
      yx = r_templates_1[0].y*signal[0].x;
      s_input_1[threadIdx.x].x = (xx - yy)*scale;
      s_input_1[threadIdx.x].y = (xy + yx)*scale;
      s_input_2[threadIdx.x].x = (xx + yy)*scale;
      s_input_2[threadIdx.x].y = (xy - yx)*scale;

      xx = r_templates_1[1].x*signal[1].x;
      yy = r_templates_1[1].y*signal[1].y;
      xy = r_templates_1[1].x*signal[1].y;
      yx = r_templates_1[1].y*signal[1].x;
      s_input_1[threadIdx.x + (KERNLEN>>2)].x = (xx - yy)*scale;
      s_input_1[threadIdx.x + (KERNLEN>>2)].y = (xy + yx)*scale;
      s_input_2[threadIdx.x + (KERNLEN>>2)].x = (xx + yy)*scale;
      s_input_2[threadIdx.x + (KERNLEN>>2)].y = (xy - yx)*scale;
		
      xx = r_templates_1[2].x*signal[2].x;
      yy = r_templates_1[2].y*signal[2].y;
      xy = r_templates_1[2].x*signal[2].y;
      yx = r_templates_1[2].y*signal[2].x;
      s_input_1[threadIdx.x + (KERNLEN>>1)].x = (xx - yy)*scale;
      s_input_1[threadIdx.x + (KERNLEN>>1)].y = (xy + yx)*scale;
      s_input_2[threadIdx.x + (KERNLEN>>1)].x = (xx + yy)*scale;
      s_input_2[threadIdx.x + (KERNLEN>>1)].y = (xy - yx)*scale;
		
      xx = r_templates_1[3].x*signal[3].x;
      yy = r_templates_1[3].y*signal[3].y;
      xy = r_templates_1[3].x*signal[3].y;
      yx = r_templates_1[3].y*signal[3].x;
      s_input_1[threadIdx.x + 3*(KERNLEN>>2)].x = (xx - yy)*scale;
      s_input_1[threadIdx.x + 3*(KERNLEN>>2)].y = (xy + yx)*scale;
      s_input_2[threadIdx.x + 3*(KERNLEN>>2)].x = (xx + yy)*scale;
      s_input_2[threadIdx.x + 3*(KERNLEN>>2)].y = (xy - yx)*scale;
		
      __syncthreads();
		
      //----------> IFFT
      do_IFFT_mk11_4elem_2vertical_no_reorder(s_input_1, s_input_2);
      //----------<
		
      // Saving data
      pos = t*useful_part_size*nConvolutions + blockIdx.x*useful_part_size + threadIdx.x;
      if( threadIdx.x>=offset && threadIdx.x<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos - offset] = pwcalc(s_input_1[threadIdx.x]);
      }
      if( (threadIdx.x+(KERNLEN>>2))>=offset && (threadIdx.x+(KERNLEN>>2))<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos + (KERNLEN>>2) - offset] = pwcalc(s_input_1[threadIdx.x + (KERNLEN>>2)]);
      }
      if( (threadIdx.x+(KERNLEN>>1))>=offset && (threadIdx.x+(KERNLEN>>1))<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos + (KERNLEN>>1) - offset] = pwcalc(s_input_1[threadIdx.x + (KERNLEN>>1)]);
      }
      if( (threadIdx.x+3*(KERNLEN>>2))>=offset && (threadIdx.x+3*(KERNLEN>>2))<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos + 3*(KERNLEN>>2) - offset] = pwcalc(s_input_1[threadIdx.x + 3*(KERNLEN>>2)]);
      }
		
		
      pos = (ZMAX-t)*useful_part_size*nConvolutions + blockIdx.x*useful_part_size + threadIdx.x;
      if( threadIdx.x>=offset && threadIdx.x<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos - offset] = pwcalc(s_input_2[threadIdx.x]);
      }
      if( (threadIdx.x+(KERNLEN>>2))>=offset && (threadIdx.x+(KERNLEN>>2))<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos + (KERNLEN>>2) - offset] = pwcalc(s_input_2[threadIdx.x + (KERNLEN>>2)]);
      }
      if( (threadIdx.x+(KERNLEN>>1))>=offset && (threadIdx.x+(KERNLEN>>1))<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos + (KERNLEN>>1) - offset] = pwcalc(s_input_2[threadIdx.x + (KERNLEN>>1)]);
      }
      if( (threadIdx.x+3*(KERNLEN>>2))>=offset && (threadIdx.x+3*(KERNLEN>>2))<(useful_part_size+offset) ) {
	d_output_plane_reduced[pos + 3*(KERNLEN>>2) - offset] = pwcalc(s_input_2[threadIdx.x + 3*(KERNLEN>>2)]);
      }
    }
	
    // now do z=0
    pos = (ZMAX>>1)*KERNLEN + threadIdx.x;
    r_templates_1[0]=__ldg(&d_templates[pos]);
    r_templates_1[1]=__ldg(&d_templates[pos + (KERNLEN>>2)]);
    r_templates_1[2]=__ldg(&d_templates[pos + (KERNLEN>>1)]);
    r_templates_1[3]=__ldg(&d_templates[pos + 3*(KERNLEN>>2)]);
	
    s_input_1[threadIdx.x].x                  = (r_templates_1[0].x*signal[0].x - r_templates_1[0].y*signal[0].y)*scale;
    s_input_1[threadIdx.x].y                  = (r_templates_1[0].x*signal[0].y + r_templates_1[0].y*signal[0].x)*scale;
    s_input_1[threadIdx.x + (KERNLEN>>2)].x   = (r_templates_1[1].x*signal[1].x - r_templates_1[1].y*signal[1].y)*scale;
    s_input_1[threadIdx.x + (KERNLEN>>2)].y   = (r_templates_1[1].x*signal[1].y + r_templates_1[1].y*signal[1].x)*scale;
    s_input_1[threadIdx.x + (KERNLEN>>1)].x   = (r_templates_1[2].x*signal[2].x - r_templates_1[2].y*signal[2].y)*scale;
    s_input_1[threadIdx.x + (KERNLEN>>1)].y   = (r_templates_1[2].x*signal[2].y + r_templates_1[2].y*signal[2].x)*scale;
    s_input_1[threadIdx.x + 3*(KERNLEN>>2)].x = (r_templates_1[3].x*signal[3].x - r_templates_1[3].y*signal[3].y)*scale;
    s_input_1[threadIdx.x + 3*(KERNLEN>>2)].y = (r_templates_1[3].x*signal[3].y + r_templates_1[3].y*signal[3].x)*scale;
	
    __syncthreads();
	
    //call inverse fft
    do_IFFT_mk11_4elem_no_reorder(s_input_1);
	
    pos = (ZMAX/2)*useful_part_size*nConvolutions + blockIdx.x*useful_part_size + threadIdx.x;
    if( threadIdx.x>=offset && threadIdx.x<(useful_part_size+offset) ) {
      d_output_plane_reduced[pos - offset] = pwcalc(s_input_1[threadIdx.x]);
    }
    if( (threadIdx.x+(KERNLEN>>2))>=offset && (threadIdx.x+(KERNLEN>>2))<(useful_part_size+offset) ) {
      d_output_plane_reduced[pos + (KERNLEN>>2) - offset] = pwcalc(s_input_1[threadIdx.x + (KERNLEN>>2)]);
    }
    if( (threadIdx.x+(KERNLEN>>1))>=offset && (threadIdx.x+(KERNLEN>>1))<(useful_part_size+offset) ) {
      d_output_plane_reduced[pos + (KERNLEN>>1) - offset] = pwcalc(s_input_1[threadIdx.x + (KERNLEN>>1)]);
    }
    if( (threadIdx.x+3*(KERNLEN>>2))>=offset && (threadIdx.x+3*(KERNLEN>>2))<(useful_part_size+offset) ) {
      d_output_plane_reduced[pos + 3*(KERNLEN>>2) - offset] = pwcalc(s_input_1[threadIdx.x + 3*(KERNLEN>>2)]);
    }
	
  }

  /** \brief Kernel wrapper function for GPU_CONV_kFFT_mk11_4elem_2v kernel function. */
  void call_kernel_GPU_CONV_kFFT_mk11_4elem_2v(const dim3 &grid_size, const dim3 &block_size, float2 const*const d_input_signal, float *const d_output_plane_reduced, float2 const*const d_templates, const int &useful_part_size, const int &offset, const int &nConvolutions, const float &scale) {
    GPU_CONV_kFFT_mk11_4elem_2v<<<grid_size, block_size>>>(d_input_signal, d_output_plane_reduced, d_templates, useful_part_size, offset, nConvolutions, scale);
  }

} //namespace astroaccelerate
  
#endif
