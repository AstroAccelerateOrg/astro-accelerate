#include "aa_device_load_data.hpp"
#include "aa_device_dedispersion_kernel.hpp"
#include "aa_device_SPS_inplace_kernel.hpp"

#include "aa_params.hpp"

namespace astroaccelerate {

  /**
   * \brief Function to load data from host memory into GPU memory.
   * \warning If the file extension of this file is *.cpp, then the code will compile but there will be a runtime CUDA error when copying to device memory.
   */
	void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short const*const host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts, float *d_dm_shifts, int nbits) {
		if(i == -1) {
			const long int length = ( t_processed + maxshift );
			const size_t size = (size_t)nchans * (size_t)length * (size_t)sizeof(unsigned short);
			//checkCudaErrors(cudaGetLastError());
			cudaMemcpy(device_pointer, host_pointer, size, cudaMemcpyHostToDevice);
			//checkCudaErrors(cudaGetLastError());
			if( (nchans>8192) && (nbits != 4) ){
				cudaMemcpy(d_dm_shifts, dmshifts, nchans*sizeof(float), cudaMemcpyHostToDevice);
				set_device_constants_dedispersion_kernel(nchans, length, t_processed);
			}
			else if ( (nbits == 4) && (nchans > 4096) ){
				cudaMemcpy(d_dm_shifts, dmshifts, nchans*sizeof(float), cudaMemcpyHostToDevice);
				set_device_constants_dedispersion_kernel(nchans, length, t_processed);
			}
			else {
				set_device_constants_dedispersion_kernel(nchans, length, t_processed, dmshifts);
			}
		}
		else if(i > 0) {
			const long int length = ( t_processed + maxshift );
			set_device_constants_dedispersion_kernel(length, t_processed);
		}
	}
   
  
	void load_chunk_data(unsigned short *device_pointer, unsigned short const*const host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts, float *d_dm_shifts, int nbits) {
		const long int length = ( t_processed + maxshift );
		const size_t size = (size_t)nchans * (size_t)length * (size_t)sizeof(unsigned short);
		//checkCudaErrors(cudaGetLastError());
		cudaMemcpy(device_pointer, host_pointer, size, cudaMemcpyHostToDevice);
		//checkCudaErrors(cudaGetLastError());
		if( (nchans>8192) && (nbits != 4) ){
			cudaMemcpy(d_dm_shifts, dmshifts, nchans*sizeof(float), cudaMemcpyHostToDevice);
			set_device_constants_dedispersion_kernel(nchans, length, t_processed);
		}
		else if ( (nbits == 4) && (nchans > 4096) ){
			cudaMemcpy(d_dm_shifts, dmshifts, nchans*sizeof(float), cudaMemcpyHostToDevice);
			set_device_constants_dedispersion_kernel(nchans, length, t_processed);
		}
		else {
			set_device_constants_dedispersion_kernel(nchans, length, t_processed, dmshifts);
		}
    }
	
    void set_dedispersion_constants(int t_processed, int maxshift){
		const long int length = ( t_processed + maxshift );
		set_device_constants_dedispersion_kernel(length, t_processed);
    }
} //namespace astroaccelerate
