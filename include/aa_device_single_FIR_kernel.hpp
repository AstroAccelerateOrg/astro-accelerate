#ifndef ASTRO_ACCELERATE_AA_DEVICE_SINGLE_FIR_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SINGLE_FIR_KERNEL_HPP

namespace astroaccelerate {

  /** \brief Kernel wrapper function for PD_FIR_GPU kernel function. */
  void call_kernel_PD_FIR_GPU(const dim3 &grid_size, const dim3 &block_size, const int &SM_size, float const *const d_input, float *const d_output, const int &nTaps, const int &nLoops, const int &nTimesamples);

  /** \brief Kernel wrapper function for PD_FIR_GPUv1 kernel function. */
  void call_kernel_PD_FIR_GPUv1(const dim3 &grid_size, const dim3 &block_size, const int &SM_size, float const *const d_input, float *const d_output, const int &nTaps, const int &nLoops, const unsigned int &nTimesamples);

  /** \brief Kernel wrapper function for Fir_L1 kernel function. */
  void call_kernel_Fir_L1(const dim3 &grid_size, const dim3 &block_size, float const *const d_input, float *const d_output, const int &nTaps, const int &nTimesamples);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SINGLE_FIR_KERNEL_HPP
