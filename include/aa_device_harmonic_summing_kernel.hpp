#ifndef ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_KERNEL_HPP

namespace astroaccelerate {
  class HRMS_ConstParams {
  public:
    static const int warp = 32;
  };
  
  class HRMS_normal : public HRMS_ConstParams {
  public:
    static const bool remove_scalloping_loss = false;
  };
  
  class HRMS_remove_scalloping_loss : public HRMS_ConstParams {
  public:
    static const bool remove_scalloping_loss = true;
  };

  /** \brief Kernel wrapper function for PHS_GPU_kernel kernel function. */
  void call_kernel_simple_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nSpectra,
      const int &nHarmonics
  );
  
  void call_kernel_greedy_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nSpectra,
      const int &nHarmonics,
      bool enable_scalloping_loss_removal
  );
  
  void call_kernel_presto_plus_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nSpectra,
      const int &nHarmonics,
      bool enable_scalloping_loss_removal
  );

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_KERNEL_HPP
