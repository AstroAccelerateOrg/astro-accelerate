#ifndef CONVOLUTION_KERNEL_H_
#define CONVOLUTION_KERNEL_H_

#include "aa_params.hpp"
#include "aa_device_convolution_kernel.hpp"
#include "aa_device_cuda_deprecated_wrappers.cuh"

namespace astroaccelerate {

	class FFT_Params {
	public:
		static const int fft_exp = -1;
		static const int fft_length = -1;
		static const int warp = 32;
	};

	class FFT_256 : public FFT_Params {
		public:
		static const int fft_exp = 8;
		static const int fft_length = 256;
		static const int fft_length_quarter = 64;
		static const int fft_length_half = 128;
		static const int fft_length_three_quarters = 192;
	};

	class FFT_512 : public FFT_Params {
		public:
		static const int fft_exp = 9;
		static const int fft_length = 512;
		static const int fft_length_quarter = 128;
		static const int fft_length_half = 256;
		static const int fft_length_three_quarters = 384;
	};

	class FFT_1024 : public FFT_Params {
		public:
		static const int fft_exp = 10;
		static const int fft_length = 1024;
		static const int fft_length_quarter = 256;
		static const int fft_length_half = 512;
		static const int fft_length_three_quarters = 768;
	};

	class FFT_2048 : public FFT_Params {
		public:
		static const int fft_exp = 11;
		static const int fft_length = 2048;
		static const int fft_length_quarter = 512;
		static const int fft_length_half = 1024;
		static const int fft_length_three_quarters = 1536;
	};

	class FFT_4096 : public FFT_Params {
		public:
		static const int fft_exp = 12;
		static const int fft_length = 4096;
		static const int fft_length_quarter = 1024;
		static const int fft_length_half = 2048;
		static const int fft_length_three_quarters = 3072;
	};



	class FFT_ConstDirection {
	public:
		static const int fft_direction = -1;
	};

	class FFT_forward : public FFT_ConstDirection {
	public:
		static const int fft_direction = 0;
	};

	class FFT_inverse : public FFT_ConstDirection {
	public:
		static const int fft_direction = 1;
	};

//-------------------------------------


	__device__ __inline__ float2 Get_W_value(int N, int m){
		float2 ctemp;
		sincosf ( -6.283185308f*fdividef( (float) m, (float) N), &ctemp.y, &ctemp.x);
		return(ctemp);
	}

	__device__ __inline__ float2 Get_W_value_inverse(int N, int m){
		float2 ctemp;
		sincosf ( 6.283185308f*fdividef( (float) m, (float) N), &ctemp.y, &ctemp.x);
		return(ctemp);
	}



	/* This function calculate the inverse FFT. It expect out of order input, that is order of the elements in the input must in the form produced by the Cooley-Tukey FFT algorithm.*/
	template<class const_params>
	__inline__ __device__ void CT_DIT_FFT_4way(float2 *s_input){
		float2 A_DFT_value, B_DFT_value, C_DFT_value, D_DFT_value;
		float2 W;
		float2 Aftemp, Bftemp, Cftemp, Dftemp;

		int local_id, warp_id;
		int j, m_param;
		int parity, itemp;
		int A_read_index, B_read_index, C_read_index, D_read_index;
		int PoT, PoTp1, q;
		
		local_id = threadIdx.x & (const_params::warp - 1);
		warp_id = threadIdx.x/const_params::warp;

		#ifdef TESTING
		int A_load_id, B_load_id, i, A_n, B_n;
		A_load_id = threadIdx.x;
		B_load_id = threadIdx.x + const_params::fft_length_quarter;
		A_n=threadIdx.x;
		B_n=threadIdx.x + const_params::fft_length_quarter;
		for(i=1; i<const_params::fft_exp; i++) {
			A_n >>= 1;
			B_n >>= 1;
			A_load_id <<= 1;
			A_load_id |= A_n & 1;
			B_load_id <<= 1;
			B_load_id |= B_n & 1;
		}
		A_load_id &= const_params::fft_length-1;
		B_load_id &= const_params::fft_length-1;
		
		//-----> Scrambling input
		A_DFT_value=s_input[A_load_id];
		B_DFT_value=s_input[A_load_id + 1];
		C_DFT_value=s_input[B_load_id];
		D_DFT_value=s_input[B_load_id + 1];
		__syncthreads();
		s_input[threadIdx.x]         = A_DFT_value;
		s_input[threadIdx.x + const_params::fft_length_half]   = B_DFT_value;
		s_input[threadIdx.x + const_params::fft_length_quarter]   = C_DFT_value;
		s_input[threadIdx.x + const_params::fft_length_three_quarters] = D_DFT_value;
		__syncthreads();
		#endif
		
		
		//-----> FFT
		//-->
		PoT=1;
		PoTp1=2;	

		//--> First iteration
		itemp=local_id&1;
		parity=(1-itemp*2);
		A_DFT_value=s_input[local_id + (warp_id<<2)*const_params::warp];
		B_DFT_value=s_input[local_id + (warp_id<<2)*const_params::warp + const_params::warp];
		C_DFT_value=s_input[local_id + (warp_id<<2)*const_params::warp + 2*const_params::warp];
		D_DFT_value=s_input[local_id + (warp_id<<2)*const_params::warp + 3*const_params::warp];
		
		__syncthreads();
		
		A_DFT_value.x=parity*A_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK, A_DFT_value.x, 1);
		A_DFT_value.y=parity*A_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK, A_DFT_value.y, 1);
		B_DFT_value.x=parity*B_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK, B_DFT_value.x, 1);
		B_DFT_value.y=parity*B_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK, B_DFT_value.y, 1);
		C_DFT_value.x=parity*C_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK, C_DFT_value.x, 1);
		C_DFT_value.y=parity*C_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK, C_DFT_value.y, 1);
		D_DFT_value.x=parity*D_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK, D_DFT_value.x, 1);
		D_DFT_value.y=parity*D_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK, D_DFT_value.y, 1);
		
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
			
			A_DFT_value.x = Aftemp.x + parity*aa_shfl_xor(AA_ASSUME_MASK, Aftemp.x, PoT);
			A_DFT_value.y = Aftemp.y + parity*aa_shfl_xor(AA_ASSUME_MASK, Aftemp.y, PoT);
			B_DFT_value.x = Bftemp.x + parity*aa_shfl_xor(AA_ASSUME_MASK, Bftemp.x, PoT);
			B_DFT_value.y = Bftemp.y + parity*aa_shfl_xor(AA_ASSUME_MASK, Bftemp.y, PoT);
			C_DFT_value.x = Cftemp.x + parity*aa_shfl_xor(AA_ASSUME_MASK, Cftemp.x, PoT);
			C_DFT_value.y = Cftemp.y + parity*aa_shfl_xor(AA_ASSUME_MASK, Cftemp.y, PoT);
			D_DFT_value.x = Dftemp.x + parity*aa_shfl_xor(AA_ASSUME_MASK, Dftemp.x, PoT);
			D_DFT_value.y = Dftemp.y + parity*aa_shfl_xor(AA_ASSUME_MASK, Dftemp.y, PoT);	
			
			PoT=PoT<<1;
			PoTp1=PoTp1<<1;
		}
		
		itemp = local_id + (warp_id<<2)*const_params::warp;
		s_input[itemp]                        = A_DFT_value;
		s_input[itemp + const_params::warp]   = B_DFT_value;
		s_input[itemp + 2*const_params::warp] = C_DFT_value;
		s_input[itemp + 3*const_params::warp] = D_DFT_value;
		
		for(q=5;q<(const_params::fft_exp-1);q++){
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


	/* This function calculate the forward FFT, but returns scrambled output, that is order of the elements in the output is not corrected after Cooley-Tukey FFT algorithm ends.*/
	template<class const_params>
	__inline__ __device__ void CT_DIF_FFT_4way(float2 *s_input){
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
		PoTm1 = const_params::fft_length_half;
		PoT   = const_params::fft_length;
		
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
		
		for(q=(const_params::fft_exp-2);q>4;q--){
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
			
			Aftemp.x = parity*A_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK, A_DFT_value.x, PoTm1);
			Aftemp.y = parity*A_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK, A_DFT_value.y, PoTm1);
			Bftemp.x = parity*B_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK, B_DFT_value.x, PoTm1);
			Bftemp.y = parity*B_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK, B_DFT_value.y, PoTm1);
			Cftemp.x = parity*C_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK, C_DFT_value.x, PoTm1);
			Cftemp.y = parity*C_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK, C_DFT_value.y, PoTm1);
			Dftemp.x = parity*D_DFT_value.x + aa_shfl_xor(AA_ASSUME_MASK, D_DFT_value.x, PoTm1);
			Dftemp.y = parity*D_DFT_value.y + aa_shfl_xor(AA_ASSUME_MASK, D_DFT_value.y, PoTm1);
			
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
		
		#ifdef TESTING
		__syncthreads();
		int A_load_id, B_load_id, i, A_n, B_n;
		A_load_id = threadIdx.x;
		B_load_id = threadIdx.x + const_params::fft_length_quarter;
		A_n=threadIdx.x;
		B_n=threadIdx.x + const_params::fft_length_quarter;
		for(i=1; i<const_params::fft_exp; i++) {
			A_n >>= 1;
			B_n >>= 1;
			A_load_id <<= 1;
			A_load_id |= A_n & 1;
			B_load_id <<= 1;
			B_load_id |= B_n & 1;
		}
		A_load_id &= const_params::fft_length-1;
		B_load_id &= const_params::fft_length-1;
		
		//-----> Scrambling input
		A_DFT_value=s_input[A_load_id];
		B_DFT_value=s_input[A_load_id + 1];
		C_DFT_value=s_input[B_load_id];
		D_DFT_value=s_input[B_load_id + 1];
		__syncthreads();
		s_input[threadIdx.x]                                           = A_DFT_value;
		s_input[threadIdx.x + const_params::fft_length_half]           = B_DFT_value;
		s_input[threadIdx.x + const_params::fft_length_quarter]        = C_DFT_value;
		s_input[threadIdx.x + const_params::fft_length_three_quarters] = D_DFT_value;
		__syncthreads();
		#endif
	}


	/* This is a kernel which performs forward FFT without reordering on data from device memory*/
	template<class const_params>
	__global__ void k_customFFT_GPU_forward(float2 *d_input, float2* d_output) {
		__shared__ float2 s_input[const_params::fft_length];
		s_input[threadIdx.x]                                            = d_input[threadIdx.x + blockIdx.x*const_params::fft_length];
		s_input[threadIdx.x + const_params::fft_length_quarter]        = d_input[threadIdx.x + blockIdx.x*const_params::fft_length + const_params::fft_length_quarter];
		s_input[threadIdx.x + const_params::fft_length_half]           = d_input[threadIdx.x + blockIdx.x*const_params::fft_length + const_params::fft_length_half];
		s_input[threadIdx.x + const_params::fft_length_three_quarters] = d_input[threadIdx.x + blockIdx.x*const_params::fft_length + const_params::fft_length_three_quarters];
		
		__syncthreads();
		CT_DIF_FFT_4way<const_params>(s_input);
		
		__syncthreads();
		d_output[threadIdx.x + blockIdx.x*const_params::fft_length]                                            = s_input[threadIdx.x];
		d_output[threadIdx.x + blockIdx.x*const_params::fft_length + const_params::fft_length_quarter]        = s_input[threadIdx.x + const_params::fft_length_quarter];
		d_output[threadIdx.x + blockIdx.x*const_params::fft_length + const_params::fft_length_half]           = s_input[threadIdx.x + const_params::fft_length_half];
		d_output[threadIdx.x + blockIdx.x*const_params::fft_length + const_params::fft_length_three_quarters] = s_input[threadIdx.x + const_params::fft_length_three_quarters];
	}


	/* Loads segment in correct way */
	template<class const_params>
	__device__ __inline__ void prepare_signal_4elem(float2* s_signal, float2 const* __restrict__ d_input_signal, int signal_length, int useful_part_size, int offset) {
		int pos = blockIdx.x*useful_part_size;
		
		pos = blockIdx.x*useful_part_size + threadIdx.x - offset;
		s_signal[threadIdx.x].x                                           = 0;
		s_signal[threadIdx.x].y                                           = 0;
		s_signal[threadIdx.x + const_params::fft_length_quarter].x        = 0;
		s_signal[threadIdx.x + const_params::fft_length_quarter].y        = 0;
		s_signal[threadIdx.x + const_params::fft_length_half].x           = 0;
		s_signal[threadIdx.x + const_params::fft_length_half].y           = 0;
		s_signal[threadIdx.x + const_params::fft_length_three_quarters].x = 0;
		s_signal[threadIdx.x + const_params::fft_length_three_quarters].y = 0;
		
		if( pos>=0 && pos<signal_length ) 
			s_signal[threadIdx.x] = d_input_signal[pos];
		
		if( (pos + const_params::fft_length_quarter)>=0 && (pos + const_params::fft_length_quarter)<signal_length ) 
			s_signal[threadIdx.x + const_params::fft_length_quarter] = d_input_signal[pos + const_params::fft_length_quarter];
			
		if( (pos + const_params::fft_length_half)>=0 && (pos + const_params::fft_length_half)<signal_length ) 	
			s_signal[threadIdx.x + const_params::fft_length_half] = d_input_signal[pos + const_params::fft_length_half];
		
		if( (pos + const_params::fft_length_three_quarters)>=0 && (pos + const_params::fft_length_three_quarters)<signal_length ) 
			s_signal[threadIdx.x + const_params::fft_length_three_quarters] = d_input_signal[pos + const_params::fft_length_three_quarters];
	}


	/* Convolution via overlap and save using custom FFT */
	template<class const_params>
	__global__ void k_GPU_conv_OLS_via_customFFT(
				float2 const* __restrict__ d_input_signal, 
				float *d_output_plane, 
				float2 const* __restrict__ d_filters, 
				int signal_length, 
				int useful_part_size, 
				int offset, 
				int nConvolutions, 
				int nFilters,
				float scale) {
		__shared__ float2 s_input_1[const_params::fft_length];
		float2 r_filter_1[4];
		float2 signal[4];
		int pos, t;
		
		// Loading signal segment
		prepare_signal_4elem<const_params>(s_input_1, d_input_signal, signal_length, useful_part_size, offset);
		
		// Forward FFT on input signal
		CT_DIF_FFT_4way<const_params>(s_input_1);
		
		// Storing FFTed signal for reuse
		signal[0]=s_input_1[threadIdx.x];
		signal[1]=s_input_1[threadIdx.x + const_params::fft_length_quarter];
		signal[2]=s_input_1[threadIdx.x + const_params::fft_length_half];
		signal[3]=s_input_1[threadIdx.x + const_params::fft_length_three_quarters];
		
		for(t=0; t<nFilters; t++){
			// Loading filters
			pos = t*const_params::fft_length + threadIdx.x;
			r_filter_1[0]=__ldg(&d_filters[pos]);
			r_filter_1[1]=__ldg(&d_filters[pos + const_params::fft_length_quarter]);
			r_filter_1[2]=__ldg(&d_filters[pos + const_params::fft_length_half]);
			r_filter_1[3]=__ldg(&d_filters[pos + const_params::fft_length_three_quarters]);

			// Convolution (complex multiplication)
			s_input_1[threadIdx.x].x                                           = r_filter_1[0].x*signal[0].x - r_filter_1[0].y*signal[0].y;
			s_input_1[threadIdx.x].y                                           = r_filter_1[0].x*signal[0].y + r_filter_1[0].y*signal[0].x;
			s_input_1[threadIdx.x + const_params::fft_length_quarter].x        = r_filter_1[1].x*signal[1].x - r_filter_1[1].y*signal[1].y;
			s_input_1[threadIdx.x + const_params::fft_length_quarter].y        = r_filter_1[1].x*signal[1].y + r_filter_1[1].y*signal[1].x;
			s_input_1[threadIdx.x + const_params::fft_length_half].x           = r_filter_1[2].x*signal[2].x - r_filter_1[2].y*signal[2].y;
			s_input_1[threadIdx.x + const_params::fft_length_half].y           = r_filter_1[2].x*signal[2].y + r_filter_1[2].y*signal[2].x;
			s_input_1[threadIdx.x + const_params::fft_length_three_quarters].x = r_filter_1[3].x*signal[3].x - r_filter_1[3].y*signal[3].y;
			s_input_1[threadIdx.x + const_params::fft_length_three_quarters].y = r_filter_1[3].x*signal[3].y + r_filter_1[3].y*signal[3].x;
			
			__syncthreads();
			
			//----------> Inverse FFT
			CT_DIT_FFT_4way<const_params>(s_input_1);
			//----------<
			
			
			// Writing out the clean part of the segment
			pos = t*useful_part_size*nConvolutions + blockIdx.x*useful_part_size + threadIdx.x;
			if( threadIdx.x>=offset && threadIdx.x<(useful_part_size+offset) ) {
				int spos = threadIdx.x;
				d_output_plane[pos - offset] = (s_input_1[spos].x*s_input_1[spos].x + s_input_1[spos].y*s_input_1[spos].y)*scale;
			}
			if( (threadIdx.x + const_params::fft_length_quarter)>=offset && (threadIdx.x + const_params::fft_length_quarter)<(useful_part_size+offset) ) {
				int spos = threadIdx.x + const_params::fft_length_quarter;
				d_output_plane[pos + const_params::fft_length_quarter - offset] = (s_input_1[spos].x*s_input_1[spos].x + s_input_1[spos].y*s_input_1[spos].y)*scale;
			}
			if( (threadIdx.x + const_params::fft_length_half)>=offset && (threadIdx.x + const_params::fft_length_half)<(useful_part_size+offset) ) {
				int spos = threadIdx.x + const_params::fft_length_half;
				d_output_plane[pos + const_params::fft_length_half - offset] = (s_input_1[spos].x*s_input_1[spos].x + s_input_1[spos].y*s_input_1[spos].y)*scale;
			}
			if( (threadIdx.x + const_params::fft_length_three_quarters)>=offset && (threadIdx.x + const_params::fft_length_three_quarters)<(useful_part_size+offset) ) {
				int spos = threadIdx.x + const_params::fft_length_three_quarters;
				d_output_plane[pos + const_params::fft_length_three_quarters - offset] = (s_input_1[spos].x*s_input_1[spos].x + s_input_1[spos].y*s_input_1[spos].y)*scale;
			}
			
			__syncthreads();
		}
	}

	//------------------------------------------------------->

	void call_kernel_k_customFFT_GPU_forward(
			const dim3 &grid_size, 
			const dim3 &block_size, 
			float2 *d_input, 
			float2* d_output, 
			int FFT_size)
		{
			
		switch(FFT_size) {
			case 256:
				k_customFFT_GPU_forward<FFT_256><<<grid_size, block_size>>>(d_input, d_output);
				break;
				
			case 512:
				k_customFFT_GPU_forward<FFT_512><<<grid_size, block_size>>>(d_input, d_output);
				break;
			
			case 1024:
				k_customFFT_GPU_forward<FFT_1024><<<grid_size, block_size>>>(d_input, d_output);
				break;

			case 2048:
				k_customFFT_GPU_forward<FFT_2048><<<grid_size, block_size>>>(d_input, d_output);
				break;
				
			case 4096:
				k_customFFT_GPU_forward<FFT_4096><<<grid_size, block_size>>>(d_input, d_output);
				break;
			
			default : 
				break;
		}
		
	}

	void call_kernel_k_GPU_conv_OLS_via_customFFT(
			const dim3 &grid_size, 
			const dim3 &block_size, 
			float2 *d_input_signal, 
			float *d_output_plane, 
			float2 *d_filters, 
			int signal_length, 
			int useful_part_size, 
			int offset, 
			int nConvolutions, 
			int nFilters,
			float scale,
			int convolution_length)
		{
			
		switch(convolution_length) {
			case 256:
				k_GPU_conv_OLS_via_customFFT<FFT_256><<<grid_size, block_size>>>(d_input_signal, d_output_plane, d_filters, signal_length, useful_part_size, offset, nConvolutions, nFilters, scale);
				break;
				
			case 512:
				k_GPU_conv_OLS_via_customFFT<FFT_512><<<grid_size, block_size>>>(d_input_signal, d_output_plane, d_filters, signal_length, useful_part_size, offset, nConvolutions, nFilters, scale);
				break;
			
			case 1024:
				k_GPU_conv_OLS_via_customFFT<FFT_1024><<<grid_size, block_size>>>(d_input_signal, d_output_plane, d_filters, signal_length, useful_part_size, offset, nConvolutions, nFilters, scale);
				break;

			case 2048:
				k_GPU_conv_OLS_via_customFFT<FFT_2048><<<grid_size, block_size>>>(d_input_signal, d_output_plane, d_filters, signal_length, useful_part_size, offset, nConvolutions, nFilters, scale);
				break;
				
			case 4096:
				k_GPU_conv_OLS_via_customFFT<FFT_4096><<<grid_size, block_size>>>(d_input_signal, d_output_plane, d_filters, signal_length, useful_part_size, offset, nConvolutions, nFilters, scale);
				break;
			
			default : 
				break;
		}
	}

} //namespace
#endif
