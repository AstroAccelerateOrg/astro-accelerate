#ifndef ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_KERNEL_HPP

namespace astroaccelerate {

  /** \brief Kernel wrapper function for PHS_GPU_kernel_old kernel function. */
void call_kernel_PHS_GPU_kernel_old(const dim3 &grid_size, const dim3 &block_size,
					   float const *const d_input, float *const d_output_SNR, ushort *const d_output_harmonics, float *const d_MSD, const int &nTimesamples, const int &nSpectra, const int &nHarmonics);

  /** \brief Kernel wrapper function for PHS_GPU_kernel kernel function. */
void call_kernel_PHS_GPU_kernel(const dim3 &grid_size, const dim3 &block_size,
				       float const *const d_input, float *const d_output_SNR, ushort *const d_output_harmonics, float *const d_MSD, const int &nTimesamples, const int &nSpectra, const int &nHarmonics);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_KERNEL_HPP
