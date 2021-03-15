#ifndef ASTRO_ACCELERATE_AA_DEVICE_SPECTRUM_WHITENING_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SPECTRUM_WHITENING_KERNEL_HPP

namespace astroaccelerate {

	/** \brief Kernel wrapper function for spectral whitening kernel function. */
	void call_kernel_spectrum_whitening_SGP1(
		const dim3 &block_size, 
		const dim3 &grid_size, 
		const int &smem_bytes, 
		const cudaStream_t &stream, 
		float *d_input, 
		unsigned long int nSamples
	);
	
	void call_kernel_segmented_MSD(
		const  dim3 &block_size, 
		const  dim3 &grid_size, 
		const  int &smem_bytes, 
		const  cudaStream_t &stream, 
		float  *d_segmented_MSD, 
		float2 *d_input, 
		int    *d_segment_sizes,
		size_t nSamples,
		int    nSegments
	);
	
	void call_kernel_segmented_median(
		const  dim3 &block_size, 
		const  dim3 &grid_size, 
		const  int &smem_bytes, 
		const  cudaStream_t &stream, 
		float  *d_segmented_MSD, 
		float2 *d_input, 
		int    *d_segment_sizes,
		size_t nSamples,
		int    nSegments
	);
	
	void call_kernel_spectrum_whitening_SGP2(
		const  dim3 &block_size, 
		const  dim3 &grid_size, 
		const  int &smem_bytes, 
		const  cudaStream_t &stream, 
		float  *d_segmented_MSD, 
		float2 *d_input,
		int    *d_segment_sizes,
		size_t nSamples,
		int    nSegments
	);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SPECTRUM_WHITENING_KERNEL_HPP
